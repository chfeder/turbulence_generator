#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2019-2022

import numpy as np
import os
from tempfile import mkstemp
import shutil
import subprocess
import re
import time
import datetime
import scipy.interpolate
import scipy.ndimage
from scipy.stats import binned_statistic, gaussian_kde
from scipy.signal import savgol_filter
import decimal
import functools
import pint
import myconstants as myconst
import matplotlibrc
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
from astropy.io import ascii
from astropy.table import Table
from colorama import Fore, Style
try:
    from ipdb import set_trace as stop # so we can call stop() to drop a script into interactive mode
except:
    from pdb import set_trace as stop # so we can call stop() to drop a script into interactive mode

# custom print() function override to show caller (module and function) information before actual print string
import builtins as __builtin__
import inspect
def print(*args, error=False, warn=False, highlight=False, color="", **kwargs):
    curframe = inspect.currentframe()
    calframe = inspect.getouterframes(curframe, 2)
    filename = calframe[1][1] # get file name
    ind = filename.rfind('/')+1 # remove leading path from filename
    filename = filename[ind:]
    funcname = calframe[1][3] # get function name
    prefix = filename+": "+funcname+": "
    ind = prefix.find("<module>: ") # remove '<module>: ' from prefix
    if ind < 0: ind = len(prefix)
    prefix = prefix[:ind]
    __builtin__.print(prefix, end="")
    message = "" # default message (can also be 'ERROR: ' or 'WARNING: ')
    color_str = "" # default color of text
    if color == "black"  : color_str = f"{Fore.BLACK}"
    if color == "red"    : color_str = f"{Fore.RED}"
    if color == "green"  : color_str = f"{Fore.GREEN}"
    if color == "yellow" : color_str = f"{Fore.YELLOW}"
    if color == "blue"   : color_str = f"{Fore.BLUE}"
    if color == "magenta": color_str = f"{Fore.MAGENTA}"
    if color == "cyan"   : color_str = f"{Fore.CYAN}"
    if color == "white"  : color_str = f"{Fore.WHITE}"
    if highlight:
        color_str = f"{Fore.GREEN}"
        if type(highlight)==int:
            if highlight==1: color_str = f"{Fore.GREEN}"
            if highlight==2: color_str = f"{Fore.CYAN}"
            if highlight==3: color_str = f"{Fore.MAGENTA}"
    if warn:
        message = "WARNING: "
        color_str = f"{Fore.YELLOW}"
    if error:
        message = "ERROR: "
        color_str = f"{Fore.RED}"
    # print
    __builtin__.print(color_str+message, end="")
    __builtin__.print(*args, **kwargs, flush=True)
    __builtin__.print(color_str, end=f"{Style.RESET_ALL}")
    # error handling (exit or stop the code)
    if error:
        if type(error)==str:
            if error == "stop":
                __builtin__.print(f"{Fore.RED}Type 'n' to return to function call and inspect.{Style.RESET_ALL} ", flush=True)
                return stop()
        exit() # exit on error
    return # return if everything is normal

# Decorator to have a function return a quantity with given unit (uses pint)
def use_unit(unit):
    use_unit.ureg = pint.UnitRegistry()
    def decorator_use_unit(func):
        @functools.wraps(func)
        def wrapper_use_unit(*args, **kwargs):
            value = func(*args, **kwargs)
            return value * use_unit.ureg(unit)
        return wrapper_use_unit
    return decorator_use_unit

# === Decorator to print the runtime of the decorated function ===
def timer_decorator(func):
    @functools.wraps(func)
    def timer_decorator_wrapper(*args, **kwargs):
        start_time = time.perf_counter() # start time
        value = func(*args, **kwargs) # call the actual function that we are decorating with @timer
        end_time = time.perf_counter() # end time
        run_time = end_time - start_time # compute time difference
        print(f"finished {func.__name__!r} in {run_time:.4f} seconds")
        return value
    return timer_decorator_wrapper

# === Decorator to print the function signature and return value ===
def debug_decorator(func):
    @functools.wraps(func)
    def debug_decorator_wrapper(*args, **kwargs):
        args_repr = [repr(a) for a in args] # get nice string representation of args
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()] # get nice string representation of keyword args
        signature = ", ".join(args_repr + kwargs_repr) # join them, so we see how the function was called
        print(f"calling {func.__name__}({signature})")
        value = func(*args, **kwargs)
        print(f"{func.__name__!r} returned {value!r}") # print the return value of the function
        return value
    return debug_decorator_wrapper

# === timer class to time different parts of code execution (named timers) ===
class timer:
    # ============= __init__ =============
    def __init__(self, name="", quiet=True):
        self.name = name # used to label the instance of timer (if needed)
        self.quiet = quiet # suppress time starting output
        self.start_time = None
        self.stop_time = None
        self.dt = None
        self.start()
    def start(self):
        self.start_time = datetime.datetime.now()
        if not self.quiet: print("timer('"+self.name+"'): start time = "+str(self.start_time))
    def stop(self):
        self.stop_time = datetime.datetime.now()
        if not self.quiet: print("timer('"+self.name+"'): stop time = "+str(self.stop_time))
    def get_dt(self):
        # check whether stop() was called; if not, call it here
        if self.stop_time is None: self.stop()
        # compute time difference in seconds
        self.dt = self.stop_time-self.start_time
    def report(self):
        if self.dt is None: self.get_dt()
        print("timer('"+self.name+"'): start = "+str(self.start_time)+", stop = "+str(self.stop_time)+
              ", dt = "+str(self.dt), ", dt (seconds) = "+str(self.dt.total_seconds())+" s")

# === executes a shell command: input string 'cmd'
def run_shell_command(cmd, quiet=False, print_only=False, **kargs):
    if (not quiet) or print_only:
        if 'color' not in kargs.keys():
            kargs['color'] = 'magenta' # set default colour for shell command print
        print(cmd, **kargs)
    if (not print_only):
        sp_result = subprocess.run(cmd, shell=True)
        return sp_result

def check_for_overwrite(filename):
    if os.path.isfile(filename):
        inp = input("Warning: file '"+filename+"' exisits; press 'p' to overwrite...")
        if inp != 'p': exit()

# === function returning all sub-directories, including hidden dirs ===
def get_dirs(dirname='.', include_base_dir=False, strip=False, verbose=1):
    # define recursive func
    def recursive_call(dirname='.', verbose=1):
        dirs = [f.path for f in os.scandir(dirname) if f.is_dir()]
        for dirname in list(dirs):
            if verbose > 1: print("adding sub-directories: ", dirs, highlight=3)
            dirs.extend(recursive_call(dirname, verbose=verbose))
        return dirs
    # call recursive func to get all sub-dirs
    dirs = recursive_call(dirname=dirname, verbose=verbose)
    if include_base_dir: dirs = ['.'] + dirs
    dirs = [x+'/' for x in dirs] # add trailing /
    if strip: dirs = [x[2:] for x in dirs]
    return dirs

# === function to search for a string pattern at the start of a line, and to replace that line ===
def replace_line_in_file(filename, search_str, new_line):
    debug = True
    fd, tempfile = mkstemp()
    with open(tempfile, 'w') as ftemp:
        with open(filename, 'r') as f:
            found_line = False
            for line in f:
                # replace
                if line.find(search_str)==0:
                    found_line = True
                    if debug==True: print(filename+": found line   : "+line.rstrip())
                    line = new_line+"\n"
                    if debug==True: print(filename+": replaced with: "+line.rstrip())
                # add lines to temporary output file
                ftemp.write(line)
    os.remove(filename)
    shutil.move(tempfile, filename)
    os.chmod(filename, 0o644)
    os.close(fd) # close file descriptor

# generates a numpy array with n uniformly distributed random numbers based on min and max
def generate_random_uniform_numbers(n=100, min=0.0, max=1.0, seed=None):
    random.seed(seed) # set the random seed; if None, random uses the system time
    random_numbers = [random.uniform(min, max) for _ in range(n)]
    return np.array(random_numbers)

# generates a numpy array with n Gaussian distributed random numbers based on mean mu and standard devitation sigma
def generate_random_gaussian_numbers(n=100, mu=0.0, sigma=1.0, seed=None):
    random.seed(seed) # set the random seed; if None, random uses the system time
    random_numbers = [random.gauss(mu, sigma) for _ in range(n)]
    return np.array(random_numbers)

# === return mean, stddev, skewness, kurtosis of input y(x) by integrating over x ===
def get_moments(x, y, xs=None, xe=None):
    if xs is None: xs = np.nanmin(x) # start of x (lower limit of integral)
    if xe is None: xe = np.nanmax(x) # end of x (upper limit of integral)
    ind = (x >= xs) & (x <= xe) # select relevant range
    xl = x[ind]
    yl = y[ind]
    # get middle values and dx
    xmid = ( xl[:-1] + xl[1:] ) / 2.0
    ymid = ( yl[:-1] + yl[1:] ) / 2.0
    dx = xl[1:] - xl[:-1]
    ret = {"mean": np.nan, "stddev": np.nan, "skew": 0.0, "kurt": 0.0} # init return dict
    # integrate to get moments and from that compute mean, stddev, skew, kurt
    norm = np.nansum(dx)
    if norm > 0:
        ret["mean"] = np.nansum(ymid*dx) / norm
        ret["stddev"] = np.sqrt(np.nansum((ymid-ret["mean"])**2*dx) / norm)
        if ret["stddev"] > 0:
            ret["skew"] = np.nansum(((ymid-ret["mean"])/ret["stddev"])**3*dx) / norm
            ret["kurt"] = np.nansum(((ymid-ret["mean"])/ret["stddev"])**4*dx) / norm - 3.0 # excess kurtosis
    return ret

# === get the PDF of data and return centred bin values ===
def get_pdf(data, range=None, bins=200):
    pdf, x = np.histogram(data, range=None, density=True, bins=bins)
    x_out = ( x[0:-1] + x[1:len(x)] ) / 2
    return pdf, x_out

# === bin data with bin_values (same size as data) into bins (number or array of bin edges) ===
def get_binned_stats(data, bin_values, bins=None, statistic='mean', **kwargs):
    binned_stats = binned_statistic(bin_values.flatten(), data.flatten(), bins=bins, statistic=statistic)
    return binned_stats.statistic, binned_stats.bin_edges

# === get the Fourier (k-space) spectrum of data with ncmp components in axis=0 ===
# === e.g., for a 64^3 dataset and 3 vector components, data.shape must be (3, 64, 64, 64) ===
# === e.g., for a 32^2 dataset with only 1 component, data.shape must be (32, 32) ===
def get_spectrum(data_in, ncmp=1):
    data = np.copy(data_in) # flatten to strip any possible extra dimensions
    if (ncmp == 1) and (data.shape[0] > 1):
        data = np.array([data]) # add an extra (fake) index, so we can index as if there were components
    num = np.array(data[0].shape) # number of points in data
    ndim = len(num) # number of dimensions
    ks = -num//2 # k start
    ke = np.array([np.min([-ks[d], num[d]//2-1]) for d in range(ndim)]) # k end
    k = get_coords(ks, ke, num, cell_centred=False) # get k vector with k=0 in the center
    if ndim == 1: k_abs = np.abs(k) # length of k vector
    if ndim  > 1: k_abs = np.sqrt((k**2).sum(axis=0)) # length of k vector
    bins = np.arange(np.max(num)//2) - 0.5 # k bins for 1D spectrum
    data_ft = []
    for d in range(ncmp):
        data_ft.append(np.fft.fftn(data[d], norm='forward')) # FFT
        data_ft[d] = np.fft.fftshift(data_ft[d]) # shift k=0 to center
    data_ft = np.array(data_ft)
    # get total power
    power_tot = (np.abs(data_ft)**2).sum(axis=0)
    # Helmholtz decomposition
    if ncmp > 1: # there is more then 1 component
        power_lgt = np.zeros(num, dtype=complex)
        if ndim == 1: power_lgt += k*data_ft[0] # 1D case: scalar product (k is a 1D array and we only use x-component data for the longitudinal power)
        if ndim >= 2: # 2D and 3D cases: scalar product (get power along k); if ndim < ncmp (i.e., 2.5D), the z-component does not enter the scalar product
            for d in range(ndim): power_lgt += k[d]*data_ft[d].T # scalar product
        power_lgt = np.abs(power_lgt/np.maximum(k_abs,1e-99))**2
        power_trv = power_tot - power_lgt
        print("tot power = "+str(power_tot.sum()/ncmp))
        print("lgt power = "+str(power_lgt.sum()/ncmp)+", relative to tot: "+str(power_lgt.sum()/power_tot.sum()))
        print("trv power = "+str(power_trv.sum()/ncmp)+", relative to tot: "+str(power_trv.sum()/power_tot.sum()))
    # bin in k shells
    spect_tot, bins = get_binned_stats(power_tot, k_abs, bins=bins)
    bin_centers = bins[:-1]+0.5
    integral_factor = bin_centers**(ndim-1)
    if ndim > 1: integral_factor *= np.pi*2*(ndim-1)
    ret = {'k': bin_centers, 'P_tot': spect_tot*integral_factor}
    if ncmp > 1:
        spect_lgt, bins = get_binned_stats(power_lgt, k_abs, bins=bins)
        spect_trv, bins = get_binned_stats(power_trv, k_abs, bins=bins)
        ret['P_lgt'] = spect_lgt*integral_factor
        ret['P_trv'] = spect_trv*integral_factor
    # return dict
    return ret

def get_kde_sample(data, n=1000, seed=1, show=False):
    kernel = gaussian_kde(data)
    data_resampled = kernel.resample(size=n, seed=seed)
    if show:
        pdf_original, x = get_pdf(data)
        pdf_resampled, x = get_pdf(data_resampled)
        plt.plot(x, pdf_original, label='original')
        plt.plot(x, kernel(x), label='KDE')
        plt.plot(x, pdf_resampled, label='resampled')
        plt.legend()
        plt.show(block=True)
    return data_resampled

# === return log10(x), but with the sign preserved ===
def log10(x):
    ret = np.copy(x)
    ind = ret != 0
    sign = np.sign(ret)
    ret[ind] = sign[ind]*np.log10(np.abs(ret[ind]))
    return ret

# === return x rounded to nfigs significant figures ===
def round(x, nfigs=3, str_ret=False):
    ret = float(np.format_float_positional(x, precision=nfigs, unique=False, fractional=False, trim='k'))
    if str_ret: ret = str(ret)
    return ret

# round a value and its error (uncertainty) to given nfigs significant figures
def round_with_error(val, val_err, nfigs=2):
    n = int(np.log10(val_err)) # displacement from ones place
    if val_err >= 1: n += 1
    scale = 10**(nfigs - n)
    val = np.round(val * scale) / scale
    val_err = np.round(val_err * scale) / scale
    return val, val_err

# === return x in E-format ===
def eform(x, prec=10, print_leading_plus=False):
    xx = decimal.Decimal(float(x))
    tup = xx.as_tuple()
    xx = xx.quantize(decimal.Decimal("1E{0}".format(len(tup[1])+tup[2]-prec-1)), decimal.ROUND_HALF_UP)
    tup = xx.as_tuple()
    exp = xx.adjusted()
    sign = '-' if tup.sign else '+' if print_leading_plus else ''
    dec = ''.join(str(i) for i in tup[1][1:prec+1])
    if prec > 0:
        return '{sign}{int}.{dec}E{exp:+03d}'.format(sign=sign, int=tup[1][0], dec=dec, exp=exp)
    elif prec == 0:
        return '{sign}{int}E{exp:+03d}'.format(sign=sign, int=tup[1][0], exp=exp)
    else:
        return None

# === initialise some standard plotting options/settings ===
def plot_init():
    import matplotlib.pyplot as plt
    from matplotlib import rc
    plt.style.use('classic')
    rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{bm}')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.close()
    return

# plot function
def plot(x=None, y=None, xerr=None, yerr=None, type=None, xlabel='x', ylabel='y', label=None, marker=None, linestyle=None, linewidth=None,
         show=False, pause=None, xlog=False, ylog=False, xlim=None, ylim=None, save=None, legend_loc='upper left', text=None, *args, **kwargs):
    if linestyle == 'long dashed': linestyle = (0,(5,5))
    if x is not None and y is not None:
        if text is not None:
            plt.text(x, y, text, *args, **kwargs)
        else:
            if type == 'scatter':
                linestyle = "None"
                if marker is None: marker = 'o'
            plt.errorbar(x, y, xerr=xerr, yerr=yerr, marker=marker, linestyle=linestyle, linewidth=linewidth, label=label, *args, **kwargs)
    if show or save:
        ax = plt.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xlim is not None: plt.xlim(xlim)
        if ylim is not None: plt.ylim(ylim)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        legend_handles, legend_labels = ax.get_legend_handles_labels()
        if legend_labels or legend_handles: plt.legend(loc=legend_loc)
        if save:
            plt.savefig(save, bbox_inches='tight')
            print(save+' written.', color='magenta')
        if show:
            block = None
            if pause: block = False
            plt.show(block=block)
            if pause:
                plt.draw()
                plt.pause(pause)
        plt.clf() # clear figure after use
    return plt

# escape latex
def tex_escape(text):
    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless{}',
        '>': r'\textgreater{}',
    }
    regex = re.compile('|'.join(re.escape(str(key)) for key in sorted(conv.keys(), key = lambda item: - len(item))))
    return regex.sub(lambda match: conv[match.group()], text)

# === read an ASCII file ===
def read_ascii(filename, astropy_read=True, read_header=True, quiet=False, max_num_lines=1e7, *args, **kwargs):
    if not quiet: print("reading data in '"+filename+"'...")
    if astropy_read: # simple, but slow
        tab = ascii.read(filename, *args, **kwargs)
    else: # manually reading and parsing the file; much faster
        with open(filename, 'r') as f:
            if read_header: header = np.array(f.readline().split()) # read header (first line)
            err = [] # error container
            dat = np.empty((int(max_num_lines),len(header))) # init output data table
            il = 0 # index to append line to output table
            for line in f: # loop through all lines in file
                try: dat[il] = np.asarray(line.split(), dtype=float); il += 1 # fill table with floats
                except: err.append(line) # append to error container
        dat = dat[:il] # resize table to correct size
        tab = Table() # make astropy table
        for i in range(len(dat[0])):
            tab[header[i]] = dat[:,i] # insert columns
    if not quiet: print("File '"+filename+"' read; (nrow,ncol) = ({:d},{:d}).".format(len(tab), len(tab.columns)))
    return tab

# === write an ASCII file ===
def write_ascii(filename, dat, format='fixed_width', delimiter="", comment=False, quiet=False, *args, **kwargs):
    ascii.write(dat, filename, overwrite=True, format=format, delimiter=delimiter, comment=comment, *args, **kwargs)
    if not quiet: print("File '"+filename+"' written; (nrow,ncol) = ({:d},{:d}).".format(len(dat), len(dat.columns)))

# smoothing/filtering data
def smooth(x, y, window_npts=11, order=3):
    xy_filtered = savgol_filter((x, y), window_npts, order)
    return xy_filtered[0], xy_filtered[1]

# === return freefall time ===
@use_unit("s")
def tff(rho):
    return np.sqrt(3.0*np.pi/(32.0*myconst.g_n*rho))

# === return Jeans length ===
def lJ(rho, c_s):
    return np.sqrt(np.pi*c_s**2/(myconst.g_n*rho))

# === return Jeans mass ===
def MJ(rho, c_s):
    return rho * 4.0*np.pi/3.0*(lJ(rho, c_s)/2.0)**3

def sink_dens_thresh(r_sink, c_s):
    return np.pi * c_s**2 / (4.0 * myconst.g_n * r_sink**2)

# === return mass ===
def mass(rho, L, spherical=False):
    ret = 0.0
    if spherical:
        ret = 4.0*np.pi/3.0 * rho * (L/2.0)**3
    else:
        ret = rho * L**3
    return ret

# === return virial parameter (spherical uniform-density approximation) ===
def alpha_vir(rho, L, sigma_v, spherical=False):
    return 5.0 * sigma_v**2 * L / (6.0 * myconst.g_n * mass(rho, L, spherical))

# === return sigma_s(Mach, b, beta) ===
def sigma_s(Mach, b=0.4, beta=1e99):
    beta_factor = beta / (beta + 1.0)
    sigma_s = np.sqrt(np.log(1.0 + b**2*Mach**2*beta_factor))
    return sigma_s

# === return sigma from input mean and mean-square ===
def get_sigma(mean, ms):
    diff = np.array(ms - mean**2)
    ind = np.where(diff < 0)[0]
    if diff.size > 0: diff[ind] = 0.0 # in case there is numeric rounding to near zero, we return 0
    return np.sqrt(diff)

# === return density Rankine-Hugoniot shock jump condition- ===
def shock_jump_rho(Mach, gamma=5.0/3.0):
    return (gamma+1) / (gamma-1+2/Mach**2)

# === return pressure Rankine-Hugoniot shock jump condition- ===
def shock_jump_p(Mach, gamma=5.0/3.0):
    return (1-gamma+2*gamma*Mach**2) / (gamma+1)

# === return temperature Rankine-Hugoniot shock jump condition- ===
def shock_jump_T(Mach, gamma=5.0/3.0):
    return shock_jump_p(Mach,gamma) / shock_jump_rho(Mach,gamma)

# return cell-centered coordinates | . | . |
#                               xmin       xmax
# or face-centred if keyword cell_centred=False
def get_1d_coords(cmin=0, cmax=1, ndim=10, cell_centred=True):
    if cell_centred:
        d = (cmax-cmin) / np.float(ndim)
        offset = d/2
    else:
        d = (cmax-cmin) / np.float(ndim-1)
        offset = 0.0
    return np.linspace(cmin+offset, cmax-offset, num=ndim)

def get_2d_coords(cmin=[0,0], cmax=[1,1], ndim=[10,10], cell_centred=True):
    cmin = np.array(cmin)
    cmax = np.array(cmax)
    ndim = np.array(ndim)
    if cmin.ndim != 1: cmin = [cmin,cmin]
    if cmax.ndim != 1: cmax = [cmax,cmax]
    if ndim.ndim != 1: ndim = [ndim,ndim]
    c0 = get_1d_coords(cmin=cmin[0], cmax=cmax[0], ndim=ndim[0], cell_centred=cell_centred)
    c1 = get_1d_coords(cmin=cmin[1], cmax=cmax[1], ndim=ndim[1], cell_centred=cell_centred)
    return np.array(np.meshgrid(c0, c1, indexing='ij'))

def get_3d_coords(cmin=[0,0,0], cmax=[1,1,1], ndim=[10,10,10], cell_centred=True):
    cmin = np.array(cmin)
    cmax = np.array(cmax)
    ndim = np.array(ndim)
    if cmin.ndim != 1: cmin = [cmin,cmin,cmin]
    if cmax.ndim != 1: cmax = [cmax,cmax,cmax]
    if ndim.ndim != 1: ndim = [ndim,ndim,ndim]
    c0 = get_1d_coords(cmin=cmin[0], cmax=cmax[0], ndim=ndim[0], cell_centred=cell_centred)
    c1 = get_1d_coords(cmin=cmin[1], cmax=cmax[1], ndim=ndim[1], cell_centred=cell_centred)
    c2 = get_1d_coords(cmin=cmin[2], cmax=cmax[2], ndim=ndim[2], cell_centred=cell_centred)
    return np.array(np.meshgrid(c0, c1, c2, indexing='ij'))

# this function takes lists or arrays as inputs,
# determining the dimensionality of the requested coordinates from the dimensionality of the inputs;
# for example, if cmin=[0,0], cmin=[1,1], ndim=[10,10], this function returns 2D corrdinates with 10 points in x=y=[0,1]
def get_coords(cmin, cmax, ndim, cell_centred=True):
    if (type(cmin) != list) and (type(cmin) != np.ndarray): print("need list or nump array inputs", error=True)
    cmin = np.array(cmin)
    cmax = np.array(cmax)
    ndim = np.array(ndim)
    if (cmin.shape != cmax.shape) or (cmax.shape != cmax.shape):
        print("Error: cmin, cmax, ndim, all must have the same shape.", error=True)
    if ndim.shape[0] == 1: return np.array(get_1d_coords(cmin[0], cmax[0], ndim[0], cell_centred))
    if ndim.shape[0] == 2: return np.array(get_2d_coords(cmin, cmax, ndim, cell_centred))
    if ndim.shape[0] == 3: return np.array(get_3d_coords(cmin, cmax, ndim, cell_centred))

# === polytropic_eos START ===
def polytropic_eos(dens, mu=2.3):
    polytropeDens1   = 0.0
    polytropeDens2   = 2.50e-16  # Masunaga & Inutsuka (2000)
    polytropeDens3   = 3.84e-13
    polytropeDens4   = 3.84e-8
    polytropeDens5   = 3.84e-3
    polytropeGamma1  = 1.0
    polytropeGamma2  = 1.1
    polytropeGamma3  = 1.4
    polytropeGamma4  = 1.1
    polytropeGamma5  = 1.666666
    polytropeKonst   = 4.0e8
    polytropeKonst1 = polytropeKonst
    polytropeKonst2 = polytropeKonst1 * polytropeDens2**(polytropeGamma1-polytropeGamma2)
    polytropeKonst3 = polytropeKonst2 * polytropeDens3**(polytropeGamma2-polytropeGamma3)
    polytropeKonst4 = polytropeKonst3 * polytropeDens4**(polytropeGamma3-polytropeGamma4)
    polytropeKonst5 = polytropeKonst4 * polytropeDens5**(polytropeGamma4-polytropeGamma5)
    densmin = [polytropeDens1, polytropeDens2, polytropeDens3, polytropeDens4, polytropeDens5, 1e99]
    Gamma = [polytropeGamma1, polytropeGamma2, polytropeGamma3, polytropeGamma4, polytropeGamma5]
    Konst = [polytropeKonst1, polytropeKonst2, polytropeKonst3, polytropeKonst4, polytropeKonst5]
    for i in range(5):
        if ((dens >= densmin[i]) and (dens <= densmin[i+1])):
            pres = Konst[i] * dens**Gamma[i]
            cs = np.sqrt(Gamma[i]*pres/dens)
            temp = pres / (dens/mu/myconst.m_p) / myconst.k_b
    ret = {'pres': pres, 'temp' : temp , 'cs' : cs}
    return ret
# === polytropic_eos END ===

# function to plot a map (of a 2D numpy array)
def show_map(image, dims=None, vmin=None, vmax=None, log=False, colorbar=True, cmap='magma', cmap_label=None, extent=None):
    # get dims
    target_dims = np.shape(image)
    if dims is not None: target_dims = dims
    # rebin
    image_rescaled = congrid(image, (target_dims[0], target_dims[1]))
    # define vmin and vmax
    if vmin is None: vmin = np.min(image_rescaled)
    if vmax is None: vmax = np.max(image_rescaled)
    # define how to normalise colours
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if log:
        # fix negative or zero lower bound
        if vmin <= 0: vmin = np.min(image_rescaled[image_rescaled > 0])
        norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    # plot the map
    plt.imshow(image_rescaled, cmap=cmap, origin='lower', interpolation='none', norm=norm, extent=extent)
    # add colorbar
    if colorbar:
        cb = plt.colorbar(label=cmap_label, pad=0.01, aspect=25)
        if not log: cb.ax.minorticks_on()
        cb.ax.yaxis.set_offset_position('left')
    # plot map
    plt.gcf().set_dpi(200)
    plt.tight_layout(pad=0.0)
    plt.show(block=True)

# === Smoothing a 2D array with a Gaussian kernel (sigma or fwhm in pixel units) ===
def gauss_smooth(input_data, sigma=None, fwhm=None, mode='wrap'):
    if sigma is None and fwhm is None:
        print("Either sigma or fwhm must be specified for Gaussian beam smoothing.")
        stop()
    if sigma is not None and fwhm is not None:
        print("Cannot set both sigma or fwhm; specify either sigma or fwhm for Gaussian beam smoothing.")
        stop()
    # work out the input sigma for the scipy smoothing function
    sigma_in = None
    if sigma is not None:
        sigma_in = np.array(sigma)
    if fwhm is not None:
        sigma_in = np.array(fwhm) / (2.0*np.sqrt(2.0*np.log(2.0))) # convert FWHM to sigma
    # extend dimensions of input sigma if needed
    if sigma_in.ndim != input_data.ndim-1:
        tmp = np.zeros(input_data.ndim)
        tmp[:] = sigma_in
        sigma_in = tmp
    # call scipy gaussian filter
    smoothed_data = scipy.ndimage.gaussian_filter(input_data, sigma=sigma_in, order=0, mode=mode)
    return smoothed_data


# ================== similar to IDL rebin (2D only) ====================
def rebin(a, outshape):
    inshape = a.shape
    # catch dimensions error
    if (len(inshape) != 2) or (len(outshape) != 2):
        print("Error: rebin currently only works for 2D arrays.")
        return
    # determine mode (compress or expand) and catch shape error
    error = ""
    mode = ""
    if outshape[0] <= inshape[0]:
        mode = "compress"
        if (inshape[0]%outshape[0] != 0) or (inshape[1]%outshape[1] != 0):
            error = "Error: rebin only works for integer multiples of input and output dimensions."
    else:
        mode = "expand"
        if (outshape[0]%inshape[0] != 0) or (outshape[1]%inshape[1] != 0):
            error = "Error: rebin only works for integer multiples of input and output dimensions."
    if error != "":
        print(error)
        return
    # if all good up to here, actually do work
    ret = 0
    if mode == "compress":
        sh = outshape[0],inshape[0]//outshape[0],outshape[1],inshape[1]//outshape[1]
        ret = a.reshape(sh).mean(-1).mean(1)
    if mode == "expand":
        print("Error: expending the array is not implemented yet.")
        return
    return ret

# ================== similar to IDL congrid (2D only) ====================
def congrid(a, outshape, method="linear"):
    inshape = a.shape
    # catch dimensions error
    if (len(inshape) != 2) or (len(outshape) != 2):
        print("Error: congrid currently only works for 2D arrays.")
        return
    xrange = lambda x: np.linspace(0.5/x, 1.0-0.5/x, x) # make it cell-centered interpolation
    a_clean = np.copy(a)
    ind_bad = np.isnan(a)
    if np.any(ind_bad):
        print("WARNING: nan values encountered. Setting them to 0 as a workaround...")
        a_clean[ind_bad] = 0.0
    f = scipy.interpolate.interp2d(xrange(inshape[1]), xrange(inshape[0]), a_clean, kind=method)
    a_new = f(xrange(outshape[0]), xrange(outshape[1]))
    return a_new


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":
    stop()
