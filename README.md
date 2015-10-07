# Compile-time fun

### Fun 0 -  Simple FEM with constexpr

Program does simple FEM compitation during compile-time. 
See [my blog](http://unthinkablecode.blogspot.com/2015/10/fun-with-c14-constexpr-part-2-simple.html) for description.

#### How to compile:

```
$ time g++ -Werror -Wall -Wextra -std=c++1y -o SimpleFEM SimpleFEM.cpp [OPTIONS]
```

You can remove 'time' before gcc but there is no point. 99% of code in this repository
is executed during compilation.

To add one of the following options add -DXXX=YYY to the command line
where XXX is option name and YYY is value.

Available options:

* ENVIRONMENT_TEMPERATURE
* HEAT_SOURCE_DENSITY (must be negative)
* CONVECTION_COEFFICIENT
* HEAT_TRANSFER_COEFFICIENT
* CROSSECTION_AREA
* ROD_LENGTH
* ELEMENT_COUNT

Beware that with ELEMENT_COUNT greater than 100 compilation time may lasts forever.

 