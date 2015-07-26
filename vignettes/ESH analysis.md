

# Analyzing Missouri River flow data for predicting 
# Emergent Sandbar Habitat

This document outlines a time series analysis of emergent sandbar habitat 
(ESH) which attempts to identify valuable flow metrics for empirical 
predictions of ESH habitat building and loss.


```r
# install.packages('flowregime')
require(flowregime)
data(siouxcity)
flows = separate_flow(siouxcity, filter_parameter = 0.925)
```


```r
plotdata = melt(fortify(flows), id.var = "Index", variable.name = "type", value.name = "flow")
plotdata$type = factor(plotdata$type, levels = c("quickflow", "baseflow"))
ggplot(plotdata, aes(x = Index, y = flow, color = type)) + geom_line() + facet_wrap(~type, 
    ncol = 1) + xlab("") + ylab("Flow (cms)") + scale_color_manual(values = c(baseflow = "blue", 
    quickflow = "red"), guide = FALSE)
```

![plot of chunk plot-flow](figure/plot-flow-1.png) 

The ESH data is contained in the geodatabase ESH_WEST.gdb; R can access file
geodatabases directly via the `rgdal` package. 


```r
# install.packages('rgdal')
eshw = NULL
require(rgdal)
```

```
## Loading required package: rgdal
## Loading required package: sp
## rgdal: version: 0.9-3, (SVN revision 530)
##  Geospatial Data Abstraction Library extensions to R successfully loaded
##  Loaded GDAL runtime: GDAL 1.11.2, released 2015/02/10
##  Path to GDAL shared files: C:/R/R-3.2.0/library/rgdal/gdal
##  GDAL does not use iconv for recoding strings.
##  Loaded PROJ.4 runtime: Rel. 4.9.1, 04 March 2015, [PJ_VERSION: 491]
##  Path to PROJ.4 shared files: C:/R/R-3.2.0/library/rgdal/proj
##  Linking to sp version: 1.1-0
```

```r
fgdbpath = "C:/GIS workspace/EFFECTS GIS/ESH_WEST.gdb"
colnames = c("IMAGE_DATE", "TOTAL_DISCHARGE_CFS", "Class_name", "Area_ha", "ESHCLASS", 
    "a_ESH", "a_baseESH")
for (l in ogrListLayers(fgdbpath)) {
    lyratt = as.data.frame(readOGR(fgdbpath, l, verbose = FALSE, stringsAsFactors = FALSE))
    lyratt["LAYER_NAME"] = l
    lyratt["LAYER_DATE"] = mean(as.Date(lyratt[["IMAGE_DATE"]]))
    lyratt["LAYER_DISCHARGE"] = mean(lyratt[["TOTAL_DISCHARGE_CFS"]])
    eshw = rbind(eshw, lyratt[c(colnames, "LAYER_NAME", "LAYER_DATE", "LAYER_DISCHARGE")])
}

require(dplyr)
```

```
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:xts':
## 
##     first, last
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
eshw.summary = eshw %>% group_by(Index = factor(LAYER_DATE), flow = round(LAYER_DISCHARGE)) %>% 
    summarize(sum_ESH = sum(a_baseESH), min_ESH = min(a_baseESH), max_ESH = max(a_baseESH), 
        mean_ESH = mean(a_baseESH), sd_ESH = sd(a_baseESH))

xb = xlim(c(min(index(flows)), max(index(flows))))
flowplot = ggplot(fortify(flows), aes(x = Index)) + geom_line(aes(y = baseflow), 
    color = "blue") + geom_line(aes(y = baseflow + quickflow), color = "red") + 
    xb + xlab("")

eshplot = ggplot(eshw.summary, aes(x = as.Date(Index), y = sum_ESH, size = flow)) + 
    geom_point(color = "black", fill = "white", pch = 21) + scale_size_continuous(breaks = c(10000, 
    20000, 30000, 40000)) + xb + xlab("") + ylab("Total ESH area (ha)") + theme(legend.position = "bottom")
grid.arrange(flowplot, eshplot, ncol = 1)
```

![plot of chunk load-esh](figure/load-esh-1.png) 
