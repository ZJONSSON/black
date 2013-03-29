# black.py
A generic Black-Scholes-Merton python calculator (MIT licence).   The main routine (black) is a generic solver that returns the value of the missing variable (usually implied volatility, but could be price, yield, etc)
```
black(cp=None, f=None, k=None, t=None, r=None, v=None, price=None, full=False, comp=inf):
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
 ```

 For examples, to solve for the implied volatility for a particular option:
 ```
 black("c",k=100,f=101,t=1,r=0.05,price=5)
 ```
The library also has convenience functions such as `impliedBlack` that has volatility missing as an input, `blackScholes` that generates a forward by combining a dividend yield with the risk free rate and `impliedBlackScholes` which has the volatility input missing.

Finally we have `impliedPair` which relies only on the strike price, time to maturity, risk-free and the prices of a call options and a put options at the strike.  This function solves for both the risk-free forward and the risk-neutral implied-volatility for this strike (European options only).  See Wilmott for further info: [Implied Volatility for same strike price for call and put](http://wilmott.com/messageview.cfm?catid=8&threadid=79623&FTVAR_MSGDBTABLE=)
