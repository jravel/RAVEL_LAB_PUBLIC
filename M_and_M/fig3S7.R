

## All subjects at V1

m.dat <- list(y=c(102, 11),
             n=c(521, 68))
str(m.dat)

(r <- m.2prop.Fn(m.dat, nChains=nChains, nItr=nItr, thin=thin))
## $logRat.ci
##     logRat       2.5%      97.5%    p-value
## -0.1682949 -0.7693194  0.3508604  0.2983264

## $fn.pval
## [1] 0.288366

## $est.rat
##    logRat
## 0.8451046

## $ML.rat
## [1] 0.8262687

## [[5]]
##           mean       2.5%     97.5%
## p[1] 0.1972283 0.16423251 0.2333074
## p[2] 0.1719663 0.09367827 0.2689715

## $p1
## [1] 0.1617647

## $p2
## [1] 0.1957774

## $m.dat
## $m.dat$y
## [1] 102  11

## $m.dat$n
## [1] 521  68
