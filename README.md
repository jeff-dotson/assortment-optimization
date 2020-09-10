Similarity-Based Product Assortment Optimization
================

## Abstract

Optimizing product assortments is a challenging problem. For example,
with a limited number of facings for a given CPG product, what products
do you put on the shelf? How many facings of each? Existing models call
for placing the most popular items on the shelf. In this project, we try
and incorporate measures of similarity between products to optimize an
assortment that accounts for product diversity.

## Model Development

We model the impact of product assortment, price, and promotional
choices on sales using a 2-stage model. In stage 1 we model the impact
on price and promotion on market share using a modified market share
attraction model. We specify the attractiveness (A) of each product as:

  
![ A\_{it} = \\beta\_{0i} + \\beta^p\_{it}p\_{it} +
\\beta^f\_{it}f\_{it} + \\beta^d\_{it} +
\\varepsilon\_{it}](https://latex.codecogs.com/png.latex?%20A_%7Bit%7D%20%3D%20%5Cbeta_%7B0i%7D%20%2B%20%20%5Cbeta%5Ep_%7Bit%7Dp_%7Bit%7D%20%2B%20%5Cbeta%5Ef_%7Bit%7Df_%7Bit%7D%20%2B%20%5Cbeta%5Ed_%7Bit%7D%20%2B%20%5Cvarepsilon_%7Bit%7D
" A_{it} = \\beta_{0i} +  \\beta^p_{it}p_{it} + \\beta^f_{it}f_{it} + \\beta^d_{it} + \\varepsilon_{it}")  

where
![\\beta\_{0i}](https://latex.codecogs.com/png.latex?%5Cbeta_%7B0i%7D
"\\beta_{0i}") is the product specific attractiveness (SKU intercept)
and ![p\_{it}](https://latex.codecogs.com/png.latex?p_%7Bit%7D "p_{it}")
is the price for product ![i](https://latex.codecogs.com/png.latex?i
"i") in time period ![t](https://latex.codecogs.com/png.latex?t "t").
![f\_{it}](https://latex.codecogs.com/png.latex?f_%7Bit%7D "f_{it}") and
![d\_{it}](https://latex.codecogs.com/png.latex?d_%7Bit%7D "d_{it}") are
binary indicators that are equal to
![1](https://latex.codecogs.com/png.latex?1 "1") if the product was
either Featured or Displayed during time period
![t](https://latex.codecogs.com/png.latex?t "t"). For reasons that will
be discussed later, we assume that:

  
![ \\varepsilon\_{it} \\sim N\\left(0, \\Sigma\_t\\right)
](https://latex.codecogs.com/png.latex?%20%5Cvarepsilon_%7Bit%7D%20%5Csim%20N%5Cleft%280%2C%20%5CSigma_t%5Cright%29%20
" \\varepsilon_{it} \\sim N\\left(0, \\Sigma_t\\right) ")  

Given this model structure we can compute the market share for an
product ![i](https://latex.codecogs.com/png.latex?i "i") as:

  
![ \\int { ... }
](https://latex.codecogs.com/png.latex?%20%5Cint%20%7B%20...%20%7D%20
" \\int { ... } ")  

Stage 2 of our approach models total category sales as:

  
![S\_t = \\alpha\_0 + \\alpha\_1 T\_t + \\alpha\_2 G\_t + \\alpha\_3
\\sum \_{ i }^{ }{ A\_{ it} } + \\alpha\_4|
R\_t|+\\eta\_t](https://latex.codecogs.com/png.latex?S_t%20%3D%20%5Calpha_0%20%2B%20%5Calpha_1%20T_t%20%2B%20%5Calpha_2%20G_t%20%2B%20%5Calpha_3%20%5Csum%20_%7B%20i%20%7D%5E%7B%20%20%7D%7B%20A_%7B%20it%7D%20%20%7D%20%2B%20%5Calpha_4%7C%20R_t%7C%2B%5Ceta_t
"S_t = \\alpha_0 + \\alpha_1 T_t + \\alpha_2 G_t + \\alpha_3 \\sum _{ i }^{  }{ A_{ it}  } + \\alpha_4| R_t|+\\eta_t")  

![S\_t](https://latex.codecogs.com/png.latex?S_t "S_t") denotes total
category sales, ![T\_t](https://latex.codecogs.com/png.latex?T_t "T_t")
is a term that captures trend and
![G\_t](https://latex.codecogs.com/png.latex?G_t "G_t") is a measure of
category seasonality in time ![t](https://latex.codecogs.com/png.latex?t
"t"). ![\\alpha\_3](https://latex.codecogs.com/png.latex?%5Calpha_3
"\\alpha_3") captures the impact of movements in total category
attractiveness on sales and should be positively valued. For example, if
all items in the category are being promoted total category
attractiveness, and by extension sales, should increase.

![|R\_t| \\in
\\left\[0,1\\right\]](https://latex.codecogs.com/png.latex?%7CR_t%7C%20%5Cin%20%5Cleft%5B0%2C1%5Cright%5D
"|R_t| \\in \\left[0,1\\right]") is the determinant of the correlation
matrix corresponding to the covariance matrix
![|\\Sigma\_t|](https://latex.codecogs.com/png.latex?%7C%5CSigma_t%7C
"|\\Sigma_t|"). The determinant of a correlation matrix lies on the
![\\left\[0,1\\right\]](https://latex.codecogs.com/png.latex?%5Cleft%5B0%2C1%5Cright%5D
"\\left[0,1\\right]") interval where
![1](https://latex.codecogs.com/png.latex?1 "1") indicates that the
correlations are equal to ![0](https://latex.codecogs.com/png.latex?0
"0"). As such,
![\\alpha\_4](https://latex.codecogs.com/png.latex?%5Calpha_4
"\\alpha_4") captures the influence of assortment variety on category
sales. We expect the value of
![\\alpha\_4](https://latex.codecogs.com/png.latex?%5Calpha_4
"\\alpha_4") to be positive, indicating that category sales increase as
the variety of the assortment increases. Given the functional form of
our model,
![e^{\\alpha\_4}](https://latex.codecogs.com/png.latex?e%5E%7B%5Calpha_4%7D
"e^{\\alpha_4}") is a scale-free measure of assortment value and can be
interpreted as a multiplier (e.g., a value \> 1 indicates that the
assortmen is increasing in variety). This will be important as we seek
to understand heterogenity in the value of assortment through the use of
a hierarchcial model.

We formall model ![S\_t](https://latex.codecogs.com/png.latex?S_t "S_t")
using a log-log demand model of the following functional form:

  
![ S\_t = \\alpha\_0 T\_t^{\\alpha\_1}G\_t^{\\alpha\_2} \\sum \_{ i }^{
}{ A\_{ it} }^{\\alpha\_3} \\alpha\_4^{| R\_t|}e^{\\eta\_t}
](https://latex.codecogs.com/png.latex?%20S_t%20%3D%20%5Calpha_0%20T_t%5E%7B%5Calpha_1%7DG_t%5E%7B%5Calpha_2%7D%20%5Csum%20_%7B%20i%20%7D%5E%7B%20%20%7D%7B%20A_%7B%20it%7D%20%20%7D%5E%7B%5Calpha_3%7D%20%5Calpha_4%5E%7B%7C%20R_t%7C%7De%5E%7B%5Ceta_t%7D%20
" S_t = \\alpha_0 T_t^{\\alpha_1}G_t^{\\alpha_2} \\sum _{ i }^{  }{ A_{ it}  }^{\\alpha_3} \\alpha_4^{| R_t|}e^{\\eta_t} ")  

For general details on GitHub usage, project organization, and project
workflow, see [Research Assistant
Training](https://github.com/marcdotson/ra-training).
