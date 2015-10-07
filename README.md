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

 | Option                    | Default |     
 | --- | --- | --- |       
 | ENVIRONMENT_TEMPERATURE   | 40      |    
 | HEAT_SOURCE_DENSITY       | -150    |    
 | CONVECTION_COEFFICIENT    | 10      |    
 | HEAT_TRANSFER_COEFFICIENT | 75      |    
 | CROSSECTION_AREA          | 1       |    
 | ROD_LENGTH                | 5       |    
 | ELEMENT_COUNT             | 50      |    

 Beware that with ELEMENT_COUNT greater than 100 compilation time may lasts forever.

 