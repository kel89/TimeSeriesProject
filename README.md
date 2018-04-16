# Time Series Project: Bank of Mexico Economic Indicators

In this project we are exploring various time series models for economic data from the Bank of Mexico.

### Report
The read/edit link for the overleaf report can be found [here](https://www.overleaf.com/15471400wcyrdyrcjrbx)

### Notes
- To get monthly data aggreated up to quarterly, in R we can use `aggregate(series, nfrequncy=4)` and it will take the mean of the three months in that quarter.


## Assignment
In this paper, you need to analyze several time series of your own choosing in R using the techniques taught in this course. This paper needs to be approximately 10 pages long, summarizing the most important results.

Your paper should at least consist of the following sections:

Data discussion: Discuss briefly what the time series represent and give the source of the data.

###Univariate Time Series Analysis:
- Start by making a plot of your time series, discuss stationarity, apply the steps of the Box-Jenkins approach (model specification, estimation, validation) to specify several appropriate models.
- Compare your models based on both in-sample and out-of-sample forecast criteria.


###Multivariate Time Series Analysis:
- Add another/several other time series to your paper and discuss the relationship among the time series that you would like to investigate. Note that a bivariate analysis is sufficient, but you are free to look at a relationship between more than 2 time series. Specify, estimate and validate some dynamic models
- Investigate Granger Causality
- Specify and estimate a VAR model + look at the impulse response functions
In addition to the previous sections, I also expect to see either a 

a GARCH analysis on an individual time series or
a cointegration analysis, estimation of VECM in case you work with I(1) time series.
Of course, adding both a GARCH and cointegration analysis is also an option.

 

Some additional guidelines:

Don’t only give R output in your paper but also discuss your results. Avoid including long descriptive paragraphs when discussing your results, but be concrete and to the point.
Multivariate Section: A nicely stated research question for 2 time series is valued more than a multivariate analysis with more than 2 time series but an unspecified/unclear research question.
At the end of the paper, include a short conclusion: what did you learn?
Do not forget to include the main R output in the paper (so within the 10 page limit)
Attach all your R code at the end of your paper in an Appendix (does not need to be within 10 page limit).