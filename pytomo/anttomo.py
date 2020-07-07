import numpy as np
import itertools as it  # only used in checkerboard tests
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from copy import deepcopy, copy
from inspect import getargspec
from matplotlib.colors import LinearSegmentedColormap, ColorConverter
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import gridspec
import os, sys, shutil, glob, pickle, shapefile, pyproj

####
# pysismo tomography bits cannibalized for aftan
####

########################################################################
# velocity map class for inverting a list of dispersion curves
########################################################################

class VelocityMap:
    """
    Class taking care of the inversion of velocities between
    pairs of stations, to produce a velocity map at a given
    period. The inversion procedure of Barmin et al. (2001)
    is applied.

    Attributes:
     - period      : period (s) of the velocity map
     - disp_curves : disp curves whose period's velocity is not nan
     - paths       : list of geodesic paths associated with pairs of stations
                     of dispersion curves
     - v0          : reference velocity (inverse of mean slowness, i.e.,
                     slowness implied by all observed travel-times)
     - dobs        : vector of observed data (differences observed-reference travel time)
     - Cinv        : inverse of covariance matrix of the data
     - G           : forward matrix, such that d = G.m
                     (m = parameter vector = (v0-v)/v at grid nodes)
     - density     : array of path densities at grid nodes
     - Q           : regularization matrix
     - mopt        : vector of best-fitting parameters, Ginv.C^-1.dobs
                     = best-fitting (v0-v)/v at grid nodes
     - covmopt     : covariance matrix of the best-fitting parameters, (Gt.C^-1.G + Q)^-1
     - Ginv        : inversion operator, (Gt.C^-1.G + Q)^-1.Gt = covmopt.Gt
     - R           : resolution matrix, (Gt.C^-1.G + Q)^-1.Gt.C^-1.G = Ginv.C^-1.G
     - Rradius     : array of radii of the cones that best-fit each line of the
                     resolution matrix

     Note that vectors (d, m) and matrixes (Cinv, G, Q, Ginv, R) are NOT
     numpy arrays, but numpy matrixes (vectors being n x 1 matrixes). This
     means that the product operation (*) on such objects is NOT the
     element-by-element product, but the real matrix product.
    """
    def __init__(self, dispersion_curves, period, vtype='phase', skipstations=(), skippairs=(),
                 resolution_fit='cone', min_resolution_height=0.1,
                 showplot=False, verbose=True, **kwargs):
        """
        Initializes the velocity map at period = *period*, from
        the observed velocities in *dispersion_curves*:
        - sets up the data vector, forward matrix and regularization matrix
        - performs the tomographic inversion to estimate the best-fitting
          parameters and the resolution matrix
        - estimates the characteristic spatial resolution by fitting a cone
          to each line of the resolution matrix

        Choose either phase or group velocity for tomography

        Specify stations and/or pairs to be skipped (if any), as lists, e.g.:
          skipstations = ['PORB', 'CAUB']
          skippairs = [('APOB', 'SPB'), ('ITAB', 'BAMB')];
        These options are useful (1) to test the influence of some given
        station(s) on the tomographic maps, and (2) to perform a 2-pass
        tomographic inversion, wherein pairs with a too large difference
        observed/predicted travel- time are excluded from the second pass.

        Select the type of function you want to fit to each resolution map
        with *resolution_fit*:
        - 'cone' to fit a cone, and report the cone's radius as characteristic
          resolution at each grid node in self.Rradius
        - 'gaussian' to fit a gaussian function, exp(-r/2.sigma^2), and report
          2.sigma as characteristic resolution at each grid node in self.Rradius

        Note that all resolutions in self.Rradius having a best-fitting
        cone height < *min_resolution_height* * max height will be
        discarded and set to nan.

        Append optional argument (**kwargs) to override default values:
        - minspectSNR       : min spectral SNR to retain velocity
        - lonstep           : longitude step of grid (default LONSTEP)
        - latstep           : latitude step of grid (default LATSTEP)
        - correlation_length: correlation length of the smoothing kernel:
                                S(r,r') = exp[-|r-r'|**2 / (2 * correlation_length**2)]
                              (default value CORRELATION_LENGTH)
        - alpha             : strength of the spatial smoothing term in the penalty
                              function (default ALPHA)
        - beta              : strength of the weighted norm penalization term in the
                              penalty function (default BETA)
        - lambda_           : parameter in the damping factor of the norm penalization
                              term, such that the norm is weighted by:
                                exp(- lambda_*path_density)
                              With a value of 0.15, penalization becomes strong when
                              path density < ~20
                              With a value of 0.30, penalization becomes strong when
                              path density < ~10
                              (default LAMBDA)

        @type dispersion_curves: list of L{DispersionCurve}
        @type skipstations: list of str
        @type skippairs: list of (str, str)
        """
        self.period = period
        self.vtype = vtype

        # reading inversion parameters
        minspectSNR = kwargs.get('minspectSNR', 5.0)
        minwavelengthfactor = kwargs.get('minwavelengthfactor', 2.5)
        lonstep = kwargs.get('lonstep', 0.3)
        latstep = kwargs.get('latstep', 0.3)
        correlation_length = kwargs.get('correlation_length', 50)
        alpha = kwargs.get('alpha', 200)
        beta = kwargs.get('beta', 50)
        lambda_ = kwargs.get('lambda_', 0.15)

        if verbose:
            print("Velocities selection criteria:")
            print("- rejecting velocities if SNR < {}".format(minspectSNR))
            print("- {} x {} deg grid".format(lonstep, latstep))
            s = "- correlation length of the smoothing kernel: {} km"
            print(s.format(correlation_length))
            print("- strength of the spatial smoothing term: {}".format(alpha))
            print("- strength of the norm penalization term: {}".format(beta))
            print("- weighting norm by exp(- {} * path_density)".format(lambda_))
            print()

        # skipping stations and pairs
        if skipstations:
            dispersion_curves = [c for c in dispersion_curves
                                 if not c.station1.name in skipstations and
                                 not c.station2.name in skipstations]
        if skippairs:
            skippairs = [set(pair) for pair in skippairs]
            dispersion_curves = [c for c in dispersion_curves
                                 if not {c.station1.name, c.station2.name} in skippairs]

        # updating parameters of dispersion curves
        new_qc = np.array([minspectSNR, minwavelengthfactor])
        for c in dispersion_curves:
            old_qc = np.array([c.minspectSNR, c.minwavelengthfactor])
            if np.any(abs(new_qc - old_qc)/old_qc) > 0.05 or \
			usewavelengthcutoff != c.usewavelengthcutoff:
                c.update_parameters(minspectSNR=minspectSNR,
                                    minwavelengthfactor=minwavelengthfactor)

        # valid dispersion curves (velocity != nan at period) and
        # associated interstation distances
        self.disp_curves = [c for c in dispersion_curves
                            if not np.isnan(c.filtered_vel_sdev_SNR(self.period,vtype=self.vtype)[0])]

        if not self.disp_curves:
            s = "No valid velocity at selected period ({} sec)"
            raise CannotPerformTomoInversion(s.format(period))

        dists = np.array([c.dist() for c in self.disp_curves])

        # getting (non nan) velocities and std devs at period
        vels, sigmav, _ = zip(*[c.filtered_vel_sdev_SNR(self.period,vtype=self.vtype)
                                for c in self.disp_curves])
        vels = np.array(vels)
        sigmav = np.array(sigmav)
        #sigmav_isnan = np.isnan(sigmav)  # this shouldn't be a problem b/c of interp?
        # (that is, if vel is not nan and we're using this curve, sdev will also not be nan)

########################################################################
# NOTE [skipping this part for now because the velocity step isn't constant anymore]
        # If the resolution in the velocities space is dv,
        # it means that a velocity v is actually anything between
        # v-dv/2 and v+dv/2, so the standard deviation cannot be
        # less than the standard dev of a uniform distribution of
        # width dv, which is dv / sqrt(12). Note that:
        #
        #   dv = max(dv_FTAN, dt_xc * v^2/dist),
        #
        # with dv_FTAN the intrinsic velocity discretization step
        # of the FTAN, and dt_xc the sampling interval of the
        # cross-correlation.

        #dv = np.maximum(FTAN_VELOCITIES_STEP, PERIOD_RESAMPLE * vels**2 / dists)
        #minsigmav = dv / np.sqrt(12)
        #sigmav[~sigmav_isnan] = np.maximum(sigmav[~sigmav_isnan],
        #                                   minsigmav[~sigmav_isnan])

        # where std dev cannot be estimated (std dev = nan),
        # assigning 3 times the mean std dev of the period
        # following Bensen et al. (2008)
        #sigmav[sigmav_isnan] = 3 * sigmav[~sigmav_isnan].mean()
########################################################################

        # ======================================================
        # setting up reference velocity and data vector
        # = vector of differences observed-reference travel time
        # ======================================================
        if verbose:
            print('Setting up reference velocity (v0) and data vector (dobs)')

        # reference velocity = inverse of mean slowness
        # mean slowness = slowness implied by observed travel-times
        #               = sum(observed travel-times) / sum(intersation distances)
        s = (dists / vels).sum() / dists.sum()
        self.v0 = 1.0 / s

        # data vector
        self.dobs = np.matrix(dists / vels - dists / self.v0).T

        # inverse of covariance matrix of the data
        if verbose:
            print('Setting up covariance matrix (C)')
        sigmad = sigmav * dists / vels**2
        self.Cinv = np.matrix(np.zeros((len(sigmav), len(sigmav))))
        np.fill_diagonal(self.Cinv, 1.0 / sigmad**2)

        # spatial grid for tomographic inversion (slightly enlarged to be
        # sure that no path will fall outside)
        lons1, lats1 = zip(*[c.station1.coord for c in self.disp_curves])
        lons2, lats2 = zip(*[c.station2.coord for c in self.disp_curves])
        tol = 0.5
        lonmin = np.floor(min(lons1 + lons2) - tol)
        nlon = np.ceil((max(lons1 + lons2) + tol - lonmin) / lonstep) + 1
        latmin = np.floor(min(lats1 + lats2) - tol)
        nlat = np.ceil((max(lats1 + lats2) + tol - latmin) / latstep) + 1
        self.grid = Grid(lonmin, lonstep, nlon, latmin, latstep, nlat)

        # geodesic paths associated with pairs of stations of dispersion curves
        if verbose:
            print('Calculating interstation paths')
        self.paths = []
        for curve, dist in zip(self.disp_curves, dists):
            # interpoint distance <= 1 km, and nb of points >= 100
            npts = max(np.ceil(dist) + 1, 100)
            path = geodesic(curve.station1.coord, curve.station2.coord, npts)
            self.paths.append(path)

        # ================================================
        # setting up forward matrix G, such that d = G.m
        #
        # G[i,j] = integral{w_j(r) / v0 ds} over path nb i
        # (w_j(r) = weight of node nb j on location r)
        # ================================================
        G = np.zeros((len(self.paths), self.grid.n_nodes()))
        if verbose:
            print('Setting up {} x {} forward matrix (G)'.format(*G.shape))
        for ipath, path in enumerate(self.paths):

            # for each point M along the path (1) we determine the Delaunay
            # triangle ABC that encloses M, (2) we locally define a cartesian
            # system on the plane ABC, (3) we locate M' (the projection of M
            # on the plane ABC) and (4) we attribute weights to A, B, C
            # corresponding to the three-point linear interpolation of A, B,
            # C at point M'.

            lon_M, lat_M = path[:, 0], path[:, 1]
            xyzM = geo2cartesian(lon_M, lat_M)

            # indexes, geographic coordinates and cartesian coordinates
            # (on unit sphere) of grid nodes of Delaunay triangle ABC
            # enclosing M
            iA, iB, iC = self.grid.indexes_delaunay_triangle(lon_M, lat_M)
            lonlatA, lonlatB, lonlatC = [self.grid.xy(index_) for index_ in (iA, iB, iC)]
            xyzA, xyzB, xyzC = [geo2cartesian(lon, lat)
                                for lon, lat in (lonlatA, lonlatB, lonlatC)]

            # projection of M on the plane ABC
            xyzMp = projection(xyzM, xyzA, xyzB, xyzC)

            # weights of nodes A, B, C in linear interpolation =
            # barycentric coordinates of M' in triangle ABC
            wA, wB, wC = barycentric_coords(xyzMp, xyzA, xyzB, xyzC)

            # attributing weights to grid nodes along path:
            # w[j, :] = w_j(r) = weights of node j along path
            nM = path.shape[0]
            w = np.zeros((self.grid.n_nodes(), nM))
            w[iA, range(nM)] = wA
            w[iB, range(nM)] = wB
            w[iC, range(nM)] = wC

            # ds = array of infinitesimal distances along path
            ds = udist(lons1=lon_M[:-1], lats1=lat_M[:-1],
                              lons2=lon_M[1:], lats2=lat_M[1:])

            # integrating w_j(r) / v0 along path using trapeze formula
            G[ipath, :] = np.sum(0.5 * (w[:, :-1] + w[:, 1:]) / self.v0 * ds, axis=-1)

        self.G = np.matrix(G)

        # path densities around grid's nodes
        if verbose:
            print("Calculating path densities")
        self.density = self.path_density()

        # =====================================================================
        # setting up regularization matrix Q = Ft.F + Ht.H
        #
        # F[i,j] = alpha * | 1                                         if i = j
        #                  | -S(ri,rj) / sum{S(ri,rj')} over j' != i]  if i!= j
        #
        # H[i,j] = beta * | exp[-lambda * path_density(ri)]  if i = j
        #                 | 0                                 if i!= j
        #
        # with S(.,.) the smoothing kernel and ri the locations grid nodes
        # =====================================================================

        # setting up distance matrix:
        # dists[i,j] = distance between nodes nb i and j
        dists = np.zeros((self.grid.n_nodes(), self.grid.n_nodes()))

        if verbose:
            print("Setting up {} x {} regularization matrix (Q)".format(*dists.shape))

        # indices of the upper right triangle of distance matrix
        # = (array of index #1, array of index #2)
        i_upper, j_upper = np.triu_indices_from(dists)
        lons_i, lats_i = self.grid.xy(i_upper)
        lons_j, lats_j = self.grid.xy(j_upper)
        # distance matrix (upper triangle)
        dists[i_upper, j_upper] = udist(lons1=lons_i, lats1=lats_i,
                                               lons2=lons_j, lats2=lats_j)
        # symmetrizing distance matrix (works because diagonal elts = 0)
        dists += dists.T

        # setting up smoothing kernel:
        # S[i,j] = K * exp[-|ri-rj|**2 / (2 * CORRELATION_LENGTH**2)]
        S = np.exp(- dists**2 / (2 * correlation_length**2))
        S /= S.sum(axis=-1) - np.diag(S)  # normalization of non-diagonal terms

        # setting up spatial regularization matrix F
        F = np.matrix(-S)
        F[np.diag_indices_from(F)] = 1
        F *= alpha

        # setting up regularization matrix Q
        # ... Ft.F part
        Q = F.T * F
        # ... Ht.H part
        for i, path_density in enumerate(self.density):
            Q[i, i] += beta**2 * np.exp(-2 * lambda_ * path_density)

        self.Q = Q

        # ===========================================================
        # setting up inversion operator Ginv = (Gt.C^-1.G + Q)^-1.Gt,
        # estimating model and setting up resolution matrix R =
        # Ginv.C^-1.G
        # ===========================================================

        # covariance matrix and inversion operator
        if verbose:
            print("Setting up covariance matrix of best-fitting params (covmopt)")
        self.covmopt = (self.G.T * self.Cinv * self.G + self.Q).I
        if verbose:
            print("Setting up inversion operator (Ginv)")
        self.Ginv = self.covmopt * self.G.T

        # vector of best-fitting parameters
        if verbose:
            print("Estimating best-fitting parameters (mopt)")
        self.mopt = self.Ginv * self.Cinv * self.dobs

        # resolution matrix
        if verbose:
            print("Setting up {0} x {0} resolution matrix (R)".format(self.G.shape[1]))
        self.R = self.Ginv * self.Cinv * self.G

        # ===========================================================
        # Estimating spatial resolution at each node of the grid,
        # Rradius.
        #
        # The i-th row of the resolution matrix, R[i,:], contains the
        # resolution map associated with the i-th grid noe, that is,
        # the estimated model we would get if there were only a point
        # velocity anomaly at node nb i. So a cone centered on node
        # nb i is fitted to the resolution map, and its radius gives
        # an indication of the spatial resolution at node nb i (i.e.,
        # the minimum distance at which two point anomalies can be
        # resolved)
        # ===========================================================

        if verbose:
            print("Estimation spatial resolution (Rradius)")

        self.Rradius = np.zeros(self.grid.n_nodes())
        heights = np.zeros(self.grid.n_nodes())

        for i, Ri in enumerate(np.array(self.R)):
            lon0, lat0 = self.grid.xy(i)

            # best-fitting cone at point (lon0, lat0)

            # Function returning the height of cone of radius *r0*
            # and peak *z0*, at a point located *r* km away from
            # the cone's center
            if resolution_fit.lower().strip() == 'cone':
                def cone_height(r, z0, r0):
                    """
                    Cone
                    """
                    return np.where(r < r0, z0 * (1 - r / r0), 0.0)
            elif resolution_fit.lower().strip() == 'gaussian':
                def cone_height(r, z0, r0):
                    """
                    Gaussian function
                    """
                    sigma = r0 / 2.0
                    return z0 * np.exp(- r**2 / (2 * sigma**2))
            else:
                s = "Unknown function to fit resolution: '{}'"
                raise Exception(s.format(resolution_fit))

            # distances between nodes and cone's center (lon0, lat0)
            lonnodes, latnodes = self.grid.xy_nodes()
            n = self.grid.n_nodes()
            rdata = udist(lons1=lonnodes, lats1=latnodes,
                                 lons2=n*[lon0], lats2=n*[lat0])

            # best possible resolution *rmin* = 2 * inter-node distance
            # -> estimating *rmin* along the meridian crossing the cone's
            #    center (conservative choice as it yields the largest
            #    possible value)
            d2rad = np.pi / 180.0
            rmin = 2 * d2rad * 6371.0 * max(self.grid.xstep * np.cos(lat0 * d2rad),
                                            self.grid.ystep)

            # fitting the above function to observed heights along nodes,
            # in array abs(Ri)
            popt, _ = curve_fit(f=cone_height, xdata=rdata, ydata=np.abs(Ri),
                                p0=[1, 2*rmin], maxfev=10000)
            z0, r0 = popt

            # reslution cannot be better than *rmin*
            r0 = max(rmin, r0)

            # appending spatial resolution to array
            self.Rradius[i] = r0
            heights[i] = z0

        self.Rradius[heights < heights.max() * min_resolution_height] = np.nan

        if showplot:
            # potting maps of velocity perturbation,
            # path density and resolution
            _ = self.plot()

    def __repr__(self):
        """
        E.g., "<Velocity map at period = 10 s>"
        """
        return '<Velocity map at period = {} s>'.format(self.period)

    def path_density(self, window=(0.5,0.5)):
        """
        Returns the path density, that is, on each node of the
        grid, the number of paths that cross the rectangular
        cell of size (window[0], window[1]) centered on
        the node.
        """
        # initializing path density
        density = np.zeros(self.grid.n_nodes())

        # coordinates of grid nodes and associated windows
        lons_nodes, lats_nodes = self.grid.xy_nodes()
        lons_min = np.expand_dims(lons_nodes - window[0] / 2.0, axis=-1)
        lons_max = np.expand_dims(lons_nodes + window[0] / 2.0, axis=-1)
        lats_min = np.expand_dims(lats_nodes - window[1] / 2.0, axis=-1)
        lats_max = np.expand_dims(lats_nodes + window[1] / 2.0, axis=-1)

        for path in self.paths:
            lons_path, lats_path = path[:, 0], path[:, 1]
            # are points of paths in windows?
            # 1st dim = grid nodes; 2nd dim = points along path
            points_in_windows = (lons_path >= lons_min) & (lons_path <= lons_max) & \
                                (lats_path >= lats_min) & (lats_path <= lats_max)
            density += np.any(points_in_windows, axis=-1)

        return density

    def traveltime_residuals(self, relative=False):
        """
        Returns the [relative] differences between predicted-observed
        travel times at each pair of stations:

          differences = predicted - observed travel-time,
                      = dpred - dobs,
          with dpred = G.mopt

          relative differences = (predicted - observed) / observed  travel-time
                               = (dpred - dobs) / (dobs + ref travel-time)

        @rtype: L{ndarray}
        """
        # flattening differences as 1D array
        diffs = np.array(self.G * self.mopt - self.dobs).flatten()
        if not relative:
            return diffs
        else:
            ttref = np.array([c.dist() / self.v0 for c in self.disp_curves])
            ttobs = np.array(self.dobs).flatten() + ttref  # observed travel-times
            return diffs / ttobs

    def velocity_residuals(self, relative=False):
        """
        Returns the [relative] differences between observed-predicted
        velocities (implied by travel times) at each pair of stations:

          differences = observed - predicted velocity,
                      = observed - predicted (dist / travel time),

        @rtype: L{matrix}
        """
        dists = np.array([c.dist() for c in self.disp_curves])
        ttref = np.array([c.dist() / self.v0 for c in self.disp_curves])
        ttobs = np.array(self.dobs).flatten() + ttref  # observed travel-times
        ttpred = np.array(self.G * self.mopt).flatten() + ttref  # predicted tt
        vobs = dists / ttobs  # observed velocities
        vpred = dists / ttpred  # predicted velocities
        if not relative:
            return vobs - vpred
        else:
            return (vobs - vpred) / vobs

    def model_norm(self):
        """
        Returns the norm of the model, ie the sum of velocities
        with v0 subtracted

          norm = abs(predicted velocity - v0).sum()

        @rtype: float
        """
        dists = np.array([c.dist() for c in self.disp_curves])
        ttref = np.array([c.dist() / self.v0 for c in self.disp_curves])
        ttpred = np.array(self.G * self.mopt).flatten() + ttref  # predicted tt
        vpred = dists / ttpred  # predicted velocities
        norm = abs(vpred - self.v0).sum()
        return norm

    def checkerboard_func(self, vmid, vmin, vmax, squaresize, shape='cos'):
        """
        Returns a checkerboard function, f(lons, lats), whose background
        value is *vmid*, and alternating min/max values are *vmin* and
        *vmax*. The centers of the anomalies are separated by *squaresize*
        (in km), and their shape is either 'gaussian' or 'cos'.

        @rtype: function
        """
        # converting square size from km to degrees
        d2rad = np.pi / 180.0
        midlat = 0.5 * (self.grid.ymin + self.grid.get_ymax())
        latwidth = squaresize / 6371.0 / d2rad
        lonwidth = squaresize / (6371.0 * np.cos(midlat * d2rad)) / d2rad

        # Basis function defining an anomaly of
        # unit height centered at (*lon0*, *lat0*).
        if shape.lower().strip() == 'gaussian':
            def basis_func(lons, lats, lon0, lat0):
                """
                Gausian anomaly , with sigma-parameter such that 3 sigma
                is the distance between the center and the border of
                the square, that is, half the distance between 2
                centers.
                """
                n = len(lons)
                r = udist(lons1=lons, lats1=lats, lons2=n*[lon0], lats2=n*[lat0])
                sigma = squaresize / 6.0
                return np.exp(- r**2 / (2 * sigma**2))
        elif shape.lower().strip() == 'cos':
            def basis_func(lons, lats, lon0, lat0):
                """
                Cosinus anomaly
                """
                x = (lons - lon0) / lonwidth
                y = (lats - lat0) / latwidth
                outside_square = (np.abs(x) >= 0.5) | (np.abs(y) >= 0.5)
                return np.where(outside_square, 0.0, np.cos(np.pi*x) * np.cos(np.pi*y))
        else:
            raise Exception("Unknown shape anomaly: " + shape)

        # coordinates of the center of the anomalies
        startlon = self.grid.xmin + lonwidth / 2.0
        stoplon = self.grid.get_xmax() + lonwidth
        centerlons = list(np.arange(startlon, stoplon, lonwidth))
        startlat = self.grid.ymin + latwidth / 2.0
        stoplat = self.grid.get_ymax() + latwidth
        centerlats = list(np.arange(startlat, stoplat, latwidth))
        centerlonlats = list(it.product(centerlons, centerlats))

        # factors by which multiply the basis function associated
        # with each center (to alternate lows and highs)
        polarities = [(centerlons.index(lon) + centerlats.index(lat)) % 2
                      for lon, lat in centerlonlats]
        factors = np.where(np.array(polarities) == 1, vmax - vmid, vmin - vmid)

        def func(lons, lats):
            """
            Checkboard function: sum of the basis functions along
            the centers defined above, times the high/low factor,
            plus background velocity.
            """
            lowhighs = [f * basis_func(lons, lats, lon0, lat0) for f, (lon0, lat0)
                        in zip(factors, centerlonlats)]
            return vmid + sum(lowhighs)

        return func

    def checkerboard_test(self, vmid, vmin, vmax, squaresize, **kwargs):
        """
        Generates synthetic data (travel time perturbations),
        dsynth, from a checkerboard model of velocities, and
        performs a tomographic inversion on them:

          m = (Gt.C^-1.G + Q)^-1.Gt.C^-1.dsynth
            = Ginv.C^-1.dsynth

        Returns the vector of best-fitting parameters, m.

        @rtype: L{matrix}
        """

        # checkerboard function
        f_checkerboard = self.checkerboard_func(vmid, vmin, vmax, squaresize, **kwargs)

        # setting up vector of synthetic data
        dsynth = np.zeros_like(self.dobs)
        for d, path, curve in zip(dsynth, self.paths, self.disp_curves):
            # array of infinitesimal distances along path
            lons, lats = path[:, 0], path[:, 1]
            ds = udist(lons1=lons[:-1], lats1=lats[:-1],
                              lons2=lons[1:], lats2=lats[1:])

            # velocities along path
            v = f_checkerboard(lons, lats)

            # travel time = integral[ds / v]
            t = np.sum(ds * 0.5 * (1.0 / v[:-1] + 1.0 / v[1:]))

            # synthetic data = travel time - ref travel time
            d[...] = t - curve.dist() / vmid

        # inverting synthetic data
        m = self.Ginv * self.Cinv * dsynth
        return m

    def plot(self, xsize=20, title=None, showplot=True, outfile=None, **kwargs):
        """
        Plots velocity perturbation, path density
        and spatial resolution, and returns the figure.

        Additional keyword args in *kwargs* are sent to
        self.plot_velocity(), self.plot_pathdensity()
        and self.plot_resolution(), when applicable

        @rtype: L{matplotlib.figure.Figure}
        """
        # bounding box
        bbox = self.grid.bbox()
        aspectratio = (bbox[3] - bbox[2]) / (bbox[1] - bbox[0])
        figsize = (xsize, aspectratio * xsize / 3.0 + 2)
        fig = plt.figure(figsize=figsize)

        # layout
        gs = gridspec.GridSpec(1, 3, wspace=0.0, hspace=0.0)

        # plotting velocity perturbation
        ax = fig.add_subplot(gs[0, 0])
        subkwargs = {'ax': ax, 'plot_title': False}
        # sending additional arguments (when applicable)
        subkwargs.update({k: kwargs[k] for k in getargspec(self.plot_velocity).args
                         if k in kwargs})
        self.plot_velocity(**subkwargs)

        # plotting path density
        ax = fig.add_subplot(gs[0, 1])
        subkwargs = {'ax': ax, 'plot_title': False, 'stationlabel': True}
        # sending additional arguments (when applicable)
        subkwargs.update({k: kwargs[k] for k in getargspec(self.plot_pathdensity).args
                         if k in kwargs})
        self.plot_pathdensity(**subkwargs)

        # plotting spatial resolution
        ax = fig.add_subplot(gs[0, 2])
        subkwargs = {'ax': ax, 'plot_title': False}
        # sending additional arguments (when applicable)
        subkwargs.update({k: kwargs[k] for k in getargspec(self.plot_resolution).args
                         if k in kwargs})
        self.plot_resolution(**subkwargs)

        # fig title
        if not title:
            # default title if not given
            title = u'Period = {} s, {} paths'
            title = title.format(self.period, len(self.paths))
        fig.suptitle(title, fontsize=16)

        gs.tight_layout(fig, rect=[0, 0, 1, 0.95])

        # saving figure
        if outfile:
            if os.path.exists(outfile):
                # backup
                shutil.copyfile(outfile, outfile + '~')
            fig.set_size_inches(figsize)
            fig.savefig(outfile, dpi=300)

        # showing figure
        if showplot:
            fig.show()

        return fig

    def plot_pathdensity(self, ax=None, xsize=10, plotdensity=True, plotpaths=True,
                         stationlabel=False, plot_title=True, showgrid=False,
                         highlight_residuals_gt=None):
        """
        Plots path density and/or interstation paths.

        Paths for which the residual observed/predicted travel-time
        is greater than *highlight_residuals_gt* (if defined) are
        highlighted as bold lines.
        """
        # bounding box
        bbox = self.grid.bbox()

        # creating figure if not given as input
        fig = None
        if not ax:
            aspectratio = (bbox[3] - bbox[2]) / (bbox[1] - bbox[0])
            # xzise has not effect if axes are given as input
            fig = plt.figure(figsize=(xsize, aspectratio * xsize), tight_layout=True)
            ax = fig.add_subplot(111)

        # plotting coasts and tectonic provinces
        basemap(ax=ax, labels=False, fill=not plotdensity, bbox=bbox)

        if plotdensity:
            # plotting path density
            d = self.grid.to_2D_array(self.density)
            extent = (self.grid.xmin, self.grid.get_xmax(),
                      self.grid.ymin, self.grid.get_ymax())
            m = ax.imshow(d.transpose(),
                          origin='bottom',
                          extent=extent,
                          interpolation='bicubic',
                          cmap=CMAP_DENSITY,
                          vmin=0)
            c = plt.colorbar(m, ax=ax, orientation='horizontal', pad=0.1)
            c.set_label('Path density')

        if plotpaths:
            # residuals observed/predicted travel-times
            res = self.traveltime_residuals() if highlight_residuals_gt else []

            # plotting paths
            for i, path in enumerate(self.paths):
                x, y = zip(*path)
                linestyle = {'color': 'grey', 'lw': 0.5}
                if highlight_residuals_gt and abs(float(res[i])) > highlight_residuals_gt:
                    # highlighting line as the travel-time error is > threshold
                    linestyle = {'color': 'black', 'lw': 1.5}
                ax.plot(x, y, '-', **linestyle)

        if showgrid:
            # plotting grid
            x, y = self.grid.xy_nodes()
            ax.plot(x, y, '+')

        # plotting stations
        self._plot_stations(ax, stationlabel=stationlabel)

        # formatting axes
        ax.set_xlim(bbox[:2])
        ax.set_ylim(bbox[2:])
        if plot_title:
            ax.set_title(u'Period = {} s, {} paths'.format(self.period, len(self.paths)))

        if fig:
            fig.show()

    def plot_velocity(self, ax=None, xsize=10, perturbation=False, plot_title=True,
                      vscale=None):
        """
        Plots velocity or perturbation relative to mean velocity
        (which is not necessarily the reference velocity)
        """
        # bounding box
        bbox = self.grid.bbox()

        # creating figure if not given as input
        fig = None
        if not ax:
            aspectratio = (bbox[3] - bbox[2]) / (bbox[1] - bbox[0])
            # xzise has not effect if axes are given as input
            fig = plt.figure(figsize=(xsize, aspectratio * xsize))
            ax = fig.add_subplot(111)

        # plotting coasts and tectonic provinces
        basemap(ax=ax, labels=False, fill=False, bbox=bbox)

        # plotting stations
        self._plot_stations(ax, stationlabel=False)

        # velocities on grid: m = (v0 - v) / v, so v = v0 / (1 + m)
        v = self.grid.to_2D_array(self.v0 / (1 + self.mopt))
        vmean = v.mean()
        if perturbation:
            # plotting % perturbation relative to mean velocity
            v = 100 * (v - vmean) / vmean

        if not vscale and perturbation:
            # symetric scale
            maxdv = np.abs(v).max()
            vscale = (-maxdv, maxdv)
        elif not vscale and not perturbation:
            # scale centered on mean velocity
            maxdv = np.abs(v - vmean).max()
            vscale = (vmean - maxdv, vmean + maxdv)

        # plotting spatial resolution as an overlay
        r = self.grid.to_2D_array(self.Rradius)

        extent = (self.grid.xmin, self.grid.get_xmax(),
                  self.grid.ymin, self.grid.get_ymax())
        m = ax.imshow(v.transpose(), origin='bottom', extent=extent,
                      interpolation='bicubic', cmap=CMAP_SEISMIC,
                      vmin=vscale[0], vmax=vscale[1])
        m1 = ax.imshow(r.transpose(), origin='bottom', extent=extent,
                      interpolation='bicubic',cmap=CMAP_MASK)
        c = plt.colorbar(m, ax=ax, orientation='horizontal', pad=0.1)
        c.set_label('Velocity perturbation (%)' if perturbation else 'Velocity (km/s)')

        # formatting axes
        ax.set_xlim(bbox[:2])
        ax.set_ylim(bbox[2:])
        if plot_title:
            ax.set_title(u'Period = {} s, {} paths'.format(self.period, len(self.paths)))

        if fig:
            fig.show()

    def plot_resolution(self, ax=None, xsize=10, plot_title=True):
        """
        Plots resolution map
        """
        # bounding box
        bbox = self.grid.bbox()

        # creating figure if not given as input
        fig = None
        if not ax:
            aspectratio = (bbox[3] - bbox[2]) / (bbox[1] - bbox[0])
            # xzise has not effect if axes are given as input
            fig = plt.figure(figsize=(xsize, aspectratio * xsize), tight_layout=True)
            ax = fig.add_subplot(111)

        # plotting coasts and tectonic provinces
        basemap(ax=ax, labels=False, fill=False, bbox=bbox)

        # plotting stations
        self._plot_stations(ax, stationlabel=False)

        # plotting spatial resolution
        r = self.grid.to_2D_array(self.Rradius)
        extent = (self.grid.xmin, self.grid.get_xmax(),
                  self.grid.ymin, self.grid.get_ymax())
        m = ax.imshow(r.transpose(), origin='bottom', extent=extent,
                      interpolation='bicubic',
                      cmap=CMAP_RESOLUTION)
        c = plt.colorbar(m, ax=ax, orientation='horizontal', pad=0.1)
        c.set_label('Spatial resolution (km)')

        # formatting axes
        ax.set_xlim(bbox[:2])
        ax.set_ylim(bbox[2:])
        if plot_title:
            ax.set_title(u'Period = {} s, {} paths'.format(self.period, len(self.paths)))

        if fig:
            fig.show()

    def plot_checkerboard(self, vmid, vmin, vmax, squaresize, axes=None, xsize=10,
                          **kwargs):
        """
        Plots checkboard model and reconstructed checkerboard
        """
        # checkerboard test
        m = self.checkerboard_test(vmid, vmin, vmax, squaresize, **kwargs)
        v = self.grid.to_2D_array(vmid / (1 + m))
        dv = 100 * (v - vmid) / vmid

        # bounding box
        bbox = self.grid.bbox()

        # creating figure if not given as input
        fig = None
        if not axes:
            aspectratio = (bbox[3] - bbox[2]) / (bbox[1] - bbox[0])
            # xzise has not effect if axes are given as input
            fig = plt.figure(figsize=(xsize, aspectratio * xsize), tight_layout=True)
            axes = [fig.add_subplot(121), fig.add_subplot(122)]

        ims = []

        # checkerboard model
        checkerboard_func = self.checkerboard_func(vmid, vmin, vmax, squaresize, **kwargs)
        lons, lats = self.grid.xy_nodes()
        a = self.grid.to_2D_array(checkerboard_func(lons, lats))
        extent = (self.grid.xmin, self.grid.get_xmax(),
                  self.grid.ymin, self.grid.get_ymax())
        im = axes[0].imshow(a.transpose(),
                            origin='bottom', extent=extent,
                            interpolation='bicubic',
                            vmin=vmin, vmax=vmax,
                            cmap=CMAP_SEISMIC)
        ims.append(im)

        # reconstructed checkerboard
        extent = (self.grid.xmin, self.grid.get_xmax(),
                  self.grid.ymin, self.grid.get_ymax())
        im = axes[1].imshow(dv.transpose(),
                            origin='bottom', extent=extent,
                            interpolation='bicubic',
                            vmin=-np.abs(dv).max(),
                            vmax=np.abs(dv).max(),
                            cmap=CMAP_SEISMIC)
        ims.append(im)

        for ax, im in zip(axes, ims):
            # coasts and tectonic provinces
            basemap(ax=ax, labels=False, fill=False, bbox=bbox)

            # stations
            self._plot_stations(ax, stationlabel=False)

            # color bar
            c = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.1)
            c.set_label('km/s' if ax is axes[0] else '% perturbation')

            # limits
            ax.set_xlim(bbox[:2])
            ax.set_ylim(bbox[2:])

        if fig:
            fig.show()

    def _plot_stations(self, ax, stationlabel):
        """
        Plots stations on map
        """
        # plotting stations
        xylabels = [c.station1.coord + (c.station1.name,) for c in self.disp_curves] + \
                   [c.station2.coord + (c.station2.name,) for c in self.disp_curves]
        xlist, ylist, labels = zip(*list(set(xylabels)))
        ax.plot(xlist, ylist, '^', color='k', ms=10, mfc='w', mew=1)

        if not stationlabel:
            return

        # stations label
        for x, y, label in zip(xlist, ylist, labels):
            ax.text(x, y, label, ha='center', va='bottom', fontsize=10, weight='bold')

########################################################################
# dispersion curve class, reading in info from aftan and snr files
########################################################################

class DispersionCurve:
    """
    Class holding a dispersion curve, i.e., velocity
    as a function of period
    """
    def __init__(self, station1, station2, inv,
                 disp2file,snrfile,trifiles=None,
                 minspectSNR=5.,
                 minwavelengthfactor=2.5):
        """
        Initiliazes the dispersion curve between the pair *station1*-*station2*
        by reading velocities, periods, and snr from provided files. List *trifiles*
        is all the trimester dispersion curves to use for calculating stdev.
        Selection parameters (used to select velocities that will participate
        to the tomographic inversion) are given in *minspectSNR* and *minwavelengthfactor*.

        @type station1: str
        @type station2: str
        @type inv: obspy.Inventory
        @type disp2file: str
        @type snrfile: str
        @type trifiles: list[str]
        """
        assert os.path.isfile(disp2file), 'no DISP file'
        assert os.path.isfile(snrfile), 'no snr file'
        # make station objects
        self.station1 = station_from_inv(station1, inv)
        self.station2 = station_from_inv(station2, inv)

        # read disp file (vgroup, vphase, iper, cper, phase)
        cper, iper, gvel, pvel = np.loadtxt(disp2file, usecols=(1,2,3,4), unpack=True)
        self.periods = iper
        self.vgroup = gvel
        self.vphase = pvel
        self.nomperiods = cper

        # read snrfile (sper, snr) and re-interpolate snr to curve periods
        sper,snr = np.loadtxt(snrfile, usecols=(0,1), unpack=True)
        nom2inst = interp1d(cper,iper,fill_value='extrapolate')
        ipsnr = nom2inst(sper)
        snr_interp = interp1d(ipsnr,snr,fill_value='extrapolate')
        self.SNRs = snr_interp(iper)

        # selection parameters
        self.minspectSNR = minspectSNR
        self.minwavelengthfactor = minwavelengthfactor

        # deal with trifiles: read and calculate stdev, interpolate to iper for full curve
        if len(trifiles) < 2 or trifiles == None:
            self.sdev = 0.2*np.ones(len(self.periods))  # TODO: something better than this

        else:
            for i in range(len(trifiles)):
                t_iper,t_gvel,t_pvel = np.loadtxt(trifiles[i], usecols=(2,3,4), unpack=True)
                gper = interp1d(t_iper,t_gvel,fill_value='extrapolate')  # interpolants for group&phase
                pper = interp1d(t_iper,t_pvel,fill_value='extrapolate')  # velocities with given periods
                if i == 0:
                    gtri = gper(self.periods)
                    ptri = pper(self.periods)
                else:
                    gtri = np.vstack((gtri,gper(self.periods)))
                    ptri = np.vstack((ptri,pper(self.periods)))

            self.sdev = np.std(ptri,axis=0)  # should this be group? Unclear
            self.g_tri_curves = gtri
            self.p_tri_curves = ptri

    def __repr__(self):
        return 'Dispersion curve between stations {}-{}'.format(self.station1.name,
                                                                self.station2.name)

    def update_parameters(self, minspectSNR=None,
                          minwavelengthfactor=None):
        """
        Updating one or more filtering parameter(s)
        """
        if not minspectSNR is None:
            self.minspectSNR = minspectSNR
        if not minwavelengthfactor is None:
            self.minwavelengthfactor = minwavelengthfactor

    def dist(self):
        """
        Interstation spacing (km)
        """
        return self.station1.dist(self.station2)

    def _get_cutoff_period(self):
        goodperiods = np.nan_to_num(self.dist() / (self.minwavelengthfactor * self.vgroup))
        ok = self.periods <= goodperiods
        return max(self.periods[ok])


    def filtered_vels(self,vtype='phase'):
        """
        Returns array of velocities. Periods not passing selection
        criteria are replaced with NaNs. 

        Selection criteria:
        1) period <= distance / (*minwavelengthfactor* * v)
        2) SNR >= *minspectSNR*

        @rtype: L{numpy.ndarray}
        """
        if self.SNRs is None:
            raise Exception("Spectral SNRs not defined")

        # cutoff period from minimum wavelength
        cutoff = self._get_cutoff_period()
        with np.errstate(invalid='ignore'):
            mask = self.periods <= cutoff

        # check for SNRs over threshold
        mask &= np.nan_to_num(self.SNRs) >= self.minspectSNR

        ## check to make sure vphase is always greater than vgroup
        #if vtype == 'phase':
        #    too_small = self.vphase <= self.vgroup
        #    mask[too_small] = False  # exclude vphase if less than vgroup

        # replacing velocities not passing the selection criteria with NaNs
        if vtype == 'group':
            return np.where(mask, self.vgroup, np.nan)
        elif vtype == 'phase':
            return np.where(mask, self.vphase, np.nan)

    def filtered_vel_sdev_SNR(self, period, vtype='phase'):
        """
        Returns a velocity, stdev, and SNR at a given period,
        or nan if the period does not satisfy the criteria for
        minimum SNR and/or minimum wavelength cutoff

        We use instantaneous periods here, which seems reasonable to me

        @type period: float
        @rtype: (float, float)
        """

        # check period against wavelength cutoff
        cutoff = self._get_cutoff_period()
        if period > cutoff:
            return np.nan, np.nan, np.nan

        # interpolate SNR and check against minimum
        snr_interp = interp1d(self.periods, self.SNRs, fill_value='extrapolate')
        snr_val = snr_interp(period)
        if snr_val < self.minspectSNR:
            return np.nan, np.nan, np.nan

        # if those checks pass, interpolate velocity and stdev
        if vtype == 'group': v_interp = interp1d(self.periods, self.vgroup, fill_value='extrapolate')
        if vtype == 'phase': v_interp = interp1d(self.periods, self.vphase, fill_value='extrapolate')
        v_val = v_interp(period)

        sdev_interp = interp1d(self.periods, self.sdev, fill_value='extrapolate')
        sdev_val = sdev_interp(period)

        return v_val, sdev_val, snr_val

########################################################################
# grid class for velocity map setup
########################################################################

class Grid:
    """
    Class holding a 2D regular rectangular spatial grid
    """
    def __init__(self, xmin, xstep, nx, ymin, ystep, ny):
        """
        Min coords, step size and nb of points of grid
        """
        self.xmin = xmin
        self.xstep = xstep
        self.nx = int(nx)
        self.ymin = ymin
        self.ystep = ystep
        self.ny = int(ny)

    def __repr__(self):
        s = '<2D grid: x = {}...{} by {}, y = {}...{} by {}>'
        return s.format(self.xmin, self.get_xmax(), self.xstep,
                        self.ymin, self.get_ymax(), self.ystep)

    def __eq__(self, other):
        """
        @type other: Grid
        """
        try:
            samegrids = (self.xmin == other.xmin and
                         self.xstep == other.xstep and
                         self.nx == other.nx and
                         self.ymin == other.ymin and
                         self.ystep == other.ystep and
                         self.ny == other.ny)
            return samegrids
        except:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_xmax(self):
        return self.xmin + (self.nx - 1) * self.xstep

    def get_ymax(self):
        return self.ymin + (self.ny - 1) * self.ystep

    def bbox(self):
        """
        Bounding box: (xmin, xmax, ymin, ymax)
        @rtype: (float, float, float, float)
        """
        return self.xmin, self.get_xmax(), self.ymin, self.get_ymax()

    def n_nodes(self):
        """
        Nb of nodes on grid
        """
        return self.nx * self.ny

    def ix_iy(self, index_):
        """
        Indexes along x and y-axis of node nb *index_*
        """
        ix = np.int_(np.array(index_) / self.ny)
        iy = np.mod(np.array(index_), self.ny)
        return ix, iy

    def xy(self, index_):
        """
        Coords of node nb *index_*
        """
        index_ = np.array(index_)

        if np.any((index_ < 0) | (index_ > self.n_nodes() - 1)):
            raise Exception('Index out of bounds')

        ix, iy = self.ix_iy(index_)
        return self._x(ix), self._y(iy)

    def xy_nodes(self):
        """
        Returns coords of all nodes of grid
        """
        return self.xy(np.arange(0, self.n_nodes()))

    def xarray(self):
        return np.linspace(self.xmin, self.get_xmax(), num=self.nx, endpoint=True)

    def yarray(self):
        return np.linspace(self.ymin, self.get_ymax(), num=self.ny, endpoint=True)

    def index_(self, ix, iy):
        """
        Index of node (ix, iy) in grid:
        - 0 : ix=0, iy=0
        - 1 : ix=0, iy=1
        - ...
        - ny: ix=1, iy=0
        - ...
        - nx*ny-1: ix=nx-1, iy=ny-1
        """
        ix = np.array(ix)
        iy = np.array(iy)

        if np.any((ix < 0) | (ix > self.nx - 1)):
            raise Exception('ix out of bounds')
        if np.any((iy < 0) | (iy > self.ny - 1)):
            raise Exception('iy out of bounds')

        return ix * self.ny + iy

    def indexes_delaunay_triangle(self, x, y):
        """
        Indexes of the grid's nodes defining the
        Delaunay triangle around point (x, y)
        """
        # x and y indexes of bottom left neighbour
        ix = self._xindex_left_neighbour(x)
        iy = self._yindex_bottom_neighbour(y)
        np.where(ix == self.nx - 1, ix - 1, ix)
        np.where(iy == self.ny - 1, iy - 1, iy)

        xratio = (x - self._x(ix)) / self.xstep
        yratio = (y - self._y(iy)) / self.ystep

        # returning indexes of vertices of bottom right triangle
        # or upper left triangle depending on location
        index1 = self.index_(ix, iy)
        index2 = np.where(xratio >= yratio, self.index_(ix+1, iy), self.index_(ix, iy+1))
        index3 = self.index_(ix+1, iy+1)

        return index1, index2, index3

    def geodetic_dist(self, index1, index2):
        """
        Geodetic distance between nodes nb *index1* and *index2*,
        whose coodinates (x, y) are treated as (lon, lat)
        """
        lon1, lat2 = self.xy(index1)
        lon2, lat2 = self.xy(index2)
        return udist(lons1=lon1, lats1=lat2, lons2=lon2, lats2=lat2)

    def to_2D_array(self, a):
        """
        Converts a sequence-like *a* to a 2D array b[ix, iy]
        such that i is the index of node (ix, iy)
        """
        b = np.zeros((self.nx, self.ny))
        ix, iy = self.ix_iy(range(self.n_nodes()))
        b[ix, iy] = np.array(a).flatten()
        return b

    def _x(self, ix):
        """
        Returns the abscissa of node nb *ix* on x-axis
        (ix = 0 ... nx-1)
        """
        ix = np.array(ix)
        if np.any((ix < 0) | (ix > self.nx - 1)):
            raise Exception('ix out of bounds')

        return self.xmin + ix * self.xstep

    def _y(self, iy):
        """
        Returns the ordinate of node nb *iy* on y-axis
        """
        iy = np.array(iy)
        if np.any((iy < 0) | (iy > self.ny - 1)):
            raise Exception('iy out of bounds')

        return self.ymin + iy * self.ystep

    def _xindex_left_neighbour(self, x):
        """
        Returns the index (along x-axis) of the grid nodes
        closest to (and on the left of) *x*
        (Index of 1st node = 0, index of last node = nx - 1)

        @rtype: Number
        """
        x = np.array(x)
        # checking bounds
        out_of_bounds = (x < self.xmin) | (x > self.get_xmax())
        if np.any(out_of_bounds):
            s = 'some x {} are out of bounds [{} - {}]'
            raise Exception(s.format(x[out_of_bounds], self.xmin, self.get_xmax()))

        # index of closest left node
        return np.int_((x - self.xmin) / self.xstep)

    def _yindex_bottom_neighbour(self, y):
        """
        Same as above method, along y axis

        @rtype: Number
        """
        y = np.array(y)
        # checking bounds
        out_of_bounds = (y < self.ymin) | (y > self.get_ymax())
        if np.any(out_of_bounds):
            s = 'some y {} are out of bounds [{} - {}]'
            raise Exception(s.format(y[out_of_bounds], self.ymin, self.get_ymax()))

        # index of closest bottom node
        return np.int_((y - self.ymin) / self.ystep)

########################################################################
# station info
########################################################################

def station_from_inv(stnm,inv,chn='LHZ'):
    """
    make a Station object [defined in this mod] out of obspy.inventory containing station
    """
    sel = inv.select(station=stnm,channel=chn)
    network = sel.networks[0].code
    coords = (sel.networks[0].stations[0].longitude, sel.networks[0].stations[0].latitude)
    return Station(stnm, network, chn, coords)

class Station:
    """
    Class to hold general station info: name, network, channel,
    and coordinates.
    """

    def __init__(self, name, network, channel,
                 coord=None):
        """
        @type name: str
        @type network: str
        @type channel: str
        @type coord: list of (float or None)
        """
        self.name = name
        self.network = network
        self.channel = channel  # only one channel allowed (should be BHZ)
        self.coord = coord if coord else (None, None)

    def __repr__(self):
        """
        e.g. <BL.10.NUPB>
        """
        return '<Station {0}.{1}.{2}>'.format(self.network, self.channel, self.name)

    def __str__(self):
        """
        @rtype: unicode
        """
        # General infos of station
        s = [u'Name    : {0}'.format(self.name),
             u'Network : {0}'.format(self.network),
             u'Channel : {0}'.format(self.channel),
             u'Lon, Lat: {0}, {1}'.format(*self.coord)]
        return u'\n'.join(s)

    def dist(self, other):
        """
        Geodesic distance (in km) between stations, using the
        WGS-84 ellipsoidal model of the Earth

        @type other: L{Station}
        @rtype: float
        """
        lon1, lat1 = self.coord
        lon2, lat2 = other.coord
        return udist(lons1=lon1, lats1=lat1, lons2=lon2, lats2=lat2)

    # =================
    # Boolean operators
    # =================
    BOOLATTRS = ['name', 'network', 'channel']

    def __eq__(self, other):
        """
        @type other: L{Station}
        """
        return all(getattr(self, att) == getattr(other, att) for att in self.BOOLATTRS)

    def __ne__(self, other):
        """
        @type other: L{Station}
        """
        return not self.__eq__(other)

    def __lt__(self, other):
        """
        @type other: L{Station}
        """
        return ([getattr(self, att) for att in self.BOOLATTRS] <
                [getattr(other, att) for att in self.BOOLATTRS])

    def __le__(self, other):
        """
        @type other: L{Station}
        """
        return self.__lt__(other) or self.__eq__(other)

    def __gt__(self, other):
        """
        @type other: L{Station}
        """
        return not self.__le__(other)

    def __ge__(self, other):
        """
        @type other: L{Station}
        """
        return not self.__lt__(other)

########################################################################
# utils, mostly for coordinate stuff
########################################################################

# reference elipsoid to calculate distance
wgs84 = pyproj.Geod(ellps='WGS84')

def udist(lons1, lats1, lons2, lats2):
    """
    Returns an array of geodetic distance(s) in km between
    points (lon1, lat1) and (lon2, lat2)
    """
    _, _, d = wgs84.inv(lons1=lons1, lats1=lats1, lons2=lons2, lats2=lats2)
    return np.array(d) / 1000.0


def geodesic(coord1, coord2, npts):
    """
    Returns a list of *npts* points along the geodesic between
    (and including) *coord1* and *coord2*, in an array of
    shape (*npts*, 2).
    @rtype: L{ndarray}
    """
    if npts < 2:
        raise Exception('nb of points must be at least 2')

    path = wgs84.npts(lon1=coord1[0], lat1=coord1[1],
                      lon2=coord2[0], lat2=coord2[1],
                      npts=npts-2)
    return np.array([coord1] + path + [coord2])


def geo2cartesian(lons, lats, r=1.0):
    """
    Converts geographic coordinates to cartesian coordinates
    """
    # spherical coordinates
    phi = np.array(lons) * np.pi / 180.0
    theta = np.pi / 2.0 - np.array(lats) * np.pi / 180.0
    # cartesian coordinates
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def projection(M, A, B, C):
    """
    Orthogonal projection of point(s) M on plane(s) ABC.
    Each point (M, A, B, C) should be a tuple of floats or
    a tuple of arrays, (x, y, z)
    """
    AB = vector(A, B)
    AC = vector(A, C)
    MA = vector(M, A)

    # unit vector u perpendicular to ABC (u = AB x AC / |AB x AC|)
    u = vectorial_product(AB, AC)
    norm_u = norm(u)
    u = [u[i] / norm_u for i in (0, 1, 2)]

    # (MA.u)u = MM' (with M' the projection of M on the plane)
    MA_dot_u = sum(MA[i] * u[i] for i in (0, 1, 2))
    MMp = [MA_dot_u * u[i] for i in (0, 1, 2)]
    xMp, yMp, zMp = [MMp[i] + M[i] for i in (0, 1, 2)]

    return xMp, yMp, zMp


def barycentric_coords(M, A, B, C):
    """
    Barycentric coordinates of point(s) M in triangle(s) ABC.
    Each point (M, A, B, C) should be a tuple of floats or
    a tuple of arrays, (x, y, z).
    Barycentric coordinate wrt A (resp. B, C) is the relative
    area of triangle MBC (resp. MAC, MAB).
    """
    MA = vector(M, A)
    MB = vector(M, B)
    MC = vector(M, C)

    # area of triangle = norm of vectorial product / 2
    wA = norm(vectorial_product(MB, MC)) / 2.0
    wB = norm(vectorial_product(MA, MC)) / 2.0
    wC = norm(vectorial_product(MA, MB)) / 2.0
    wtot = wA + wB + wC

    return wA / wtot, wB / wtot, wC / wtot


def vector(A, B):
    """
    Vector(s) AB. A and B should be tuple of floats or
    tuple of arrays, (x, y, z).
    """
    return tuple(np.array(B[i]) - np.array(A[i]) for i in (0, 1, 2))


def vectorial_product(u, v):
    """
    Vectorial product u x v. Vectors u, v should be tuple of
    floats or tuple of arrays, (ux, uy, uz) and (vx, vy, vz)
    """
    return (u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0])


def norm(u):
    """
    Norm of vector(s) u, which should be a tuple of
    floats or a tuple of arrays, (ux, uy, uz).
    """
    return np.sqrt(u[0]**2 + u[1]**2 + u[2]**2)

########################################################################
# plot basemap with shapefiles
########################################################################

def basemap(coastfile=None, tectfile=None, tectlabels=None, tectcolors=None, \
            ax=None, labels=True, axeslabels=True, fill=True, bbox=None):
    """
    Plots base map: coasts (file *coastfile*), tectonic provinces
    file  *tectfile*) and labels (file *tectlabels*). Labels are
    plotted if *labels* = True. Tectonic provinces are filled
    (according to colors in dict *tectcolors*) if *fill* = True.
    """
    fig = None
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    # plotting coasts
    if coastfile:
        sf = shapefile.Reader(coastfile)
        for shape in sf.shapes():
            # adding polygon(s)
            parts = list(shape.parts) + [len(shape.points)]
            partlims = zip(parts[:-1], parts[1:])
            for i1, i2 in partlims:
                points = shape.points[i1:i2]
                x, y = zip(*points)
                ax.plot(x, y, '-', lw=0.75, color='k')

    # plotting tectonic provinces
    if tectfile:
        sf = shapefile.Reader(tectfile)
        for sr in sf.shapeRecords():
            tectcategory = sr.record[0]
            color = next((tectcolors[k] for k in tectcolors.keys()
                         if k in tectcategory), 'white')
            shape = sr.shape
            parts = list(shape.parts) + [len(shape.points)]
            partlims = zip(parts[:-1], parts[1:])
            if fill:
                polygons = [Polygon(shape.points[i1:i2]) for i1, i2 in partlims]
                tectprovince = PatchCollection(polygons, facecolor=color,
                                               edgecolor='0.663', linewidths=0.5)
                ax.add_collection(tectprovince)
            else:
                for i1, i2 in partlims:
                    x, y = zip(*shape.points[i1:i2])
                    ax.plot(x, y, '-', color='0.663', lw=0.5)

    if labels and tectlabels:
        # plotting tectonic labels within bounding box
        sf = shapefile.Reader(tectlabels)
        for sr in sf.shapeRecords():
            label, angle = sr.record
            label = label.replace('\\', '\n')
            label = label.replace('Guapore', u'Guaporé').replace('Sao', u'São')
            x, y = sr.shape.points[0]
            if not bbox or bbox[0] < x < bbox[1] and bbox[2] < y < bbox[3]:
                ax.text(x, y, label, ha='center', va='center', color='grey',
                        fontsize=10, weight='bold', rotation=angle)

    # setting up axes
    ax.set_aspect('equal')
    if axeslabels:
        ax.set_xlabel('longitude (deg)')
        ax.set_ylabel('latitude (deg)')
        ax.grid(True)
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(False)
    if bbox:
        ax.set_xlim(bbox[:2])
        ax.set_ylim(bbox[2:])

    if fig:
        fig.show()


########################################################################
# Constants and parameters
########################################################################

EPS = 1.0E-6

########################################################################
# colormap things
########################################################################

# custom color map for seismic anomalies
# --------------------------------------
c = ColorConverter()
colors = ['black', 'red', 'gold', 'white',
          'white', 'aquamarine', 'blue', 'magenta']
values = [-1.0, -0.35, -0.1, -0.025,
          0.025, 0.1, 0.35, 1.0]
rgblist = [c.to_rgb(s) for s in colors]
reds, greens, blues = zip(*rgblist)
cdict = {}
for x, r, g, b in zip(values, reds, greens, blues):
    v = (x - min(values)) / (max(values) - min(values))
    cdict.setdefault('red', []).append((v, r, r))
    cdict.setdefault('green', []).append((v, g, g))
    cdict.setdefault('blue', []).append((v, b, b))
CMAP_SEISMIC = LinearSegmentedColormap('customseismic', cdict)

# custom color map for spatial resolution
# ---------------------------------------
colors = ['black', 'red', 'yellow', 'green', 'white']
values = [0, 0.25, 0.5, 0.75,  1.0]
rgblist = [c.to_rgb(s) for s in colors]
reds, greens, blues = zip(*rgblist)
cdict = {}
for x, r, g, b in zip(values, reds, greens, blues):
    v = (x - min(values)) / (max(values) - min(values))
    cdict.setdefault('red', []).append((v, r, r))
    cdict.setdefault('green', []).append((v, g, g))
    cdict.setdefault('blue', []).append((v, b, b))
CMAP_RESOLUTION = LinearSegmentedColormap('customresolution', cdict)
CMAP_RESOLUTION.set_bad(color='0.85')

# custom color map for masking resolution
# ---------------------------------------
colors = ['white','white']
values = [0, 1.0]
alphas = [0.,0.]
rgblist = [c.to_rgb(s) for s in colors]
reds, greens, blues = zip(*rgblist)
cdict = {}
for x, r, g, b, al in zip(values, reds, greens, blues, alphas):
    v = (x - min(values)) / (max(values) - min(values))
    cdict.setdefault('red', []).append((v, r, r))
    cdict.setdefault('green', []).append((v, g, g))
    cdict.setdefault('blue', []).append((v, b, b))
    cdict.setdefault('alpha',[]).append((v,al,al))
CMAP_MASK = LinearSegmentedColormap('customresolution', cdict)
CMAP_MASK.set_bad(color='k',alpha=.3)

# custom color map for path density
# ---------------------------------------
colors = ['white', 'cyan', 'green', 'yellow', 'red', 'black']
values = [0, 0.05, 0.1, 0.25, 0.5,  1.0]
rgblist = [c.to_rgb(s) for s in colors]
reds, greens, blues = zip(*rgblist)
cdict = {}
for x, r, g, b in zip(values, reds, greens, blues):
    v = (x - min(values)) / (max(values) - min(values))
    cdict.setdefault('red', []).append((v, r, r))
    cdict.setdefault('green', []).append((v, g, g))
    cdict.setdefault('blue', []).append((v, b, b))
CMAP_DENSITY = LinearSegmentedColormap('customdensity', cdict)

#########################################################################
# errors
#########################################################################

class CannotPerformTomoInversion(Exception):
    """
    Cannot perform tomographic inversion (e.g., because no data is available)
    """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
