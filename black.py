# -*- coding: iso-8859-1 -*-
# black.py (C) ziggy.jonsson.nyc@gmail.com (MIT Licence)

from __future__ import division
from scipy.stats import norm
from scipy.optimize import fsolve
from numpy import inf
from math import log, sqrt, exp

from math import *


def black(cp=None, f=None, k=None, t=None, r=None, v=None, price=None, full=False, comp=inf):
    """ General Purpose BSM machine. Leave one parameter blank and that is the parameter that will
        be solved for.
            cp   = "c" for a call, anything else assumes a put
            f    = Forward Price of the underlying asset
            k    = Strike Price
            t    = Time until maturity (in years)
            r    = Interest rate
            v    = Implied volatility
            comp = How many times interest is compounded a year (default = inf, i.e. continous rates)
            full = If True, function returns a dictionary with price and all sensitivities
                   Otherwise only the calculated/calibrated parameter will be returned
    """

    D1 = lambda f, k, sigma, t: (log(f/k)+(sigma**2/2)*t)/(sigma*sqrt(t))
    D2 = lambda f, k, sigma, t: D1(f, k, sigma, t)-sigma*sqrt(t)

    if comp != inf:
        r = log((1+r/comp)**comp)   # convert rates
    optionType = 1 if cp.upper()[0] == "C" else -1   # If the first letter is not "c" or "C" it must be a put option
    parameters = locals().copy()

    def _black(cp, f, k, t, r, v, price, full, calibration=False, **kwargs):
        d1, d2 = D1(f, k, v, t), D2(f, k, v, t)
        price = optionType * exp(-r*t)*(f*norm.cdf(optionType*d1)-k*norm.cdf(optionType*d2))

        if not calibration and full:
            return {
                "price": price,
                "delta": optionType * norm.cdf(d1),
                "gamma": norm.pdf(d1)/(f*v*sqrt(t)),
                "vega": f*norm.pdf(d1)*sqrt(t)
            }
        else:
            return price

    def _calibrate(value, field, parameters):
        parameters.update({field: value[0]})
        return abs(_black(calibration=True, **parameters) - price, )

    missing = [a for a in ["f", "k", "t", "r", "v", "price"] if parameters[a] is None]

    if len(missing) > 1:
        raise Exception("Too many missing variables from: %s " % missing)
    if len(missing) == 0:
        raise Exception("All variables assigned - nothing to solve")

    if missing[0] != "price":
        # if we are missing any parameter different from price we need to calibrate
        result = fsolve(_calibrate, 0.1, args=(missing[0], parameters))

        if full is False:
            return result   # If full=False we simply return the calibrated parameter

    return _black(**parameters)


def impliedBlack(cp, f, k, t, r, price, comp=inf):
    """ Solves for implied volatility, given all other parameters
            cp   = "c" for a call, anything else assumes a put
            f    = Forward Price of the underlying asset
            k    = Strike Price
            t    = Time until maturity (in years)
            r    = Interest rate
            price= Current Option Price
            comp = How many times interest is compounded a year (default = inf, i.e. continous rates)
            full = If True, function returns a dictionary with price and all sensitivities
                   Otherwise only the calculated/calibrated parameter will be returned
    """
    if comp != inf:
        r = log((1+r/comp)**comp)   # convert rates

    def _calibrate(v, price, cp, f, k, t, r):
        return abs(black(cp, f, k, t, r, v, comp=inf)-price)
    return fsolve(_calibrate, 0.1, args=(price, cp, f, k, t, r))


def blackScholes(cp, s, k, t, r, d, v, full=False, comp=inf):
    """ blackScholes, risk-free-rate and dividend-yield make up the forward rate
            cp   = "c" for a call, anything else assumes a put
            f    = Forward Price of the underlying asset
            k    = Strike Price
            t    = Time until maturity (in years)
            r    = Interest rate
            d    = Dividend yield
            v    = Implied volatility
            comp = How many times interest is compounded a year (default = inf, i.e. continous rates)
            full = If True, function returns a dictionary with price and all sensitivities
                   Otherwise only the calculated/calibrated parameter will be returned
    """
    if comp != inf:
        r = toExp(r, comp)
        d = toExp(d, comp)

    f = s*exp(t*(d-r))
    results = black(cp, f, k, t, r, v, full, comp=inf)

    if not full:
        return results
    else:
        for key, item in results.items():
            if key in ["delta", "gamma"]:
                item /= f/s
        results["fwd"] = f
        return results


def impliedBlackScholes(cp, s, k, t, r, d, price, comp=inf):
    """ Solves for implied volatility, given all other parameters
            cp   = "c" for a call, anything else assumes a put
            s    = Spot Price of the underlying asset
            k    = Strike Price
            t    = Time until maturity (in years)
            r    = Interest rate
            d    = Dividend yield
            v    = Implied volatility
            comp = How many times interest is compounded a year (default = inf, i.e. continous rates)
            full = If True, function returns a dictionary with price and all sensitivities
                   Otherwise only the calculated/calibrated parameter will be returned
    """
    def _calibrate(v, price, cp, f, k, t, r):
        return abs(blackScholes(cp, s, k, t, r, d, v, comp=inf)-price)

    if comp != inf:
        r = toExp(r, comp)
        d = toExp(d, comp)

    return fsolve(_calibrate, 0.1, args=(price, cp, s, k, t, r))


def impliedPair(k, t, r, callPrice, putPrice, full=False, comp=inf):
    """ Given put/call prices we solve for the forward and the implied volatility
            k    = Strike Price
            t    = Time until maturity (in years)
            r    = Interest rate
            callPrice = Price of a call @ strike
            putPrice = Price of a put @ strike
            comp = How many times interest is compounded a year (default = inf, i.e. continous rates)
            full = If True, function returns a dictionary with price and all sensitivities
                   Otherwise only the calculated/calibrated parameter will be returned
    """
    if comp != inf:
        r = log((1+r/comp)**comp)   # convert rates
    f = (callPrice-putPrice)*exp(r*t)+k
    v = impliedBlack("c", f, k, t, r, callPrice, comp=comp)
    if full:
        results = black("c", f, k, t, r, v, full=True, comp=comp)
        results["fwd"] = f
        results["vol"] = v
        return results
    else:
        return v
