---
title: CCHS Ontario population for microsimulations
author: Hana Fu
date: "21 March, 2022"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    number_sections: true
---
  







<!-- # Create factor variables -->
  

```r
cchs2013_14_Ontario$sex <- factor(cchs2013_14_Ontario$DHH_SEX, labels=c("Male","Female"))

cchs2013_14_Ontario$agecatold <- cut(cchs2013_14_Ontario$DHH_AGE, c(12,30,40,50,60,70,80,105), 
                                  include.lowest = TRUE, right = FALSE)
levels(cchs2013_14_Ontario$agecatold) <- c("12-29","30-39","40-49","50-59","60-69","70-79","80+")

cchs2013_14_Ontario$agecat <- cut(cchs2013_14_Ontario$DHH_AGE, c(12,16,25,35,55,105), 
                                  include.lowest = TRUE, right = FALSE)
levels(cchs2013_14_Ontario$agecat) <- c("12-15","16-24","25-34","35-54","55+")

cchs2013_14_Ontario$educat <- factor(cchs2013_14_Ontario$EDUDR04, 
                                     labels=c("Less than secondary","Secondary school grad", 
                                              "Some post-secondary", "Post-secondary grad", "NS"),
                                     levels =c(1,2,3,4,9))

cchs2013_14_Ontario$incpercat <- with(cchs2013_14_Ontario,
                                      ifelse(INCDPER < 10 | INCDPER %in% c(96,99), INCDPER,
                                             ifelse(INCDPER %in% c(10,11), 10,
                                                    ifelse(INCDPER %in% c(12,13), 11, 
                                                           ifelse(INCDPER == 14, 12, NA)))))

cchs2013_14_Ontario$incpercat <- factor(cchs2013_14_Ontario$incpercat,
                                        labels = c("No income", "LESS THAN $5,000", "$5,000 TO $9,999", "$10,000 TO $14,999", "$15,000 TO $19,999", "$20,000 TO $29,999", "$30,000 TO $39,999", "$40,000 TO $49,999", "$50,000 TO $59,999", "$60,000 TO $79,999","$80,000 TO 99,999","$100k or more",'NA',"NS"),
                                        levels =c(1:12,96,99))

cchs2013_14_Ontario$inchhcat <- with(cchs2013_14_Ontario,
                                     ifelse(INCDHH < 10 | INCDHH %in% c(96,99), INCDHH,
                                            ifelse(INCDHH %in% c(10,11), 10,
                                                   ifelse(INCDHH %in% c(12,13), 11, 
                                                          ifelse(INCDHH %in% c(14,15), 12, NA)))))
cchs2013_14_Ontario$inchhcat <- factor(cchs2013_14_Ontario$inchhcat,
                                       labels = c("No income", "LESS THAN $5,000", "$5,000 TO $9,999", "$10,000 TO $14,999", "$15,000 TO $19,999", "$20,000 TO $29,999", "$30,000 TO $39,999", "$40,000 TO $49,999", "$50,000 TO $59,999", "$60,000 TO $79,999","$80,000 TO 99,999","$100k or more","NS"),
                                       levels =c(1:12,99))

cchs2013_14_Ontario$ethncat <- factor(cchs2013_14_Ontario$SDCDCGT,
                                      labels = c('White','Black','Korean','Filipino','Japanese','Chinese','South Asian','South East Asian','Arab','West Asian','Latin American','Other Racial or Cultural Origin','Multiple Racial/Cultural Origins','NA','NS'),
                                      levels = c(1:13,96,99))

cchs2013_14_Ontario$maritalcat <- factor(cchs2013_14_Ontario$DHH_MS,
                                         labels = c('Married','Common-Law','Widowed','Separated','Divorced','Single/Never Married','DK','Refusal','NS'),
                                         levels = c(1:6,97,98,99))

cchs2013_14_Ontario$geourcat <- factor(cchs2013_14_Ontario$GEODUR2,
                                       labels=c("Population Centre", "Rural"),
                                       levels=c(1,2))

cchs2013_14_Ontario$alcfreqcat <- factor(cchs2013_14_Ontario$ALC_2,
                                         labels = c('less than once a month','once a month','2 to 3 times a month','once a week','2 to 3 times a week','4 to 6 times a week','every day','NA','DK','Refusal','NS'),
                                         levels = c(1:7,96:99))

cchs2013_14_Ontario$alctypecat <- factor(cchs2013_14_Ontario$ALCDTTM,
                                         labels = c("Regular drinker","Occasional drinker","Did not drink in the last 12 months","NS"),
                                         levels = c(1,2,3,9))

cchs2013_14_Ontario$alchvycat <- with(cchs2013_14_Ontario,
                                      ifelse(ALC_3 %in% 1:2, 2, # 1 is heavy drinking
                                             ifelse(ALC_3 %in% 3:6, 1, # 2 did not have heaving drinking
                                                    ifelse(ALC_3 %in% 96:99, ALC_3, NA))))

cchs2013_14_Ontario$alchvycat <- factor(cchs2013_14_Ontario$alchvycat,
                                        labels = c("Heavy drinking in the past year at least once per month","Did not heavy drink in the past year","NA","DK","Refusal","NS"),
                                        levels = c(1:2,96:99))
```

# Implement survey design


```r
# See link for instructions https://stackoverflow.com/questions/45642026/survey-weights-and-boostrap-wieghts-to-get-counts-and-cis

library(survey)
cchsDesign <- svrepdesign(data = cchs2013_14_Ontario,
                        weights = ~WTS_S,
                        repweights = "BSW[0-5]+",
                        type = "bootstrap",
                        combined.weights = TRUE,
                        mse = TRUE)
```


# Comparison of original and input data for microsimulation

CCHS baseline data with survey weights are compared with data with rounded weights used for microsimulations. Comparison is based on the proportion of people belonging to different categories of each variable of interest. The two are virtually the same, as observed in the table and barplots below. They are however usedful for determining whether some categories should be merged together due to very low proportions, such as household income categories of low income.



<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Variable </th>
   <th style="text-align:left;"> Category </th>
   <th style="text-align:right;"> CCHS weighted proportion </th>
   <th style="text-align:right;"> Rounded proportions for microsim. </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Age </td>
   <td style="text-align:left;"> 12-15 </td>
   <td style="text-align:right;"> 0.0500778 </td>
   <td style="text-align:right;"> 0.0500750 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age </td>
   <td style="text-align:left;"> 16-24 </td>
   <td style="text-align:right;"> 0.1437500 </td>
   <td style="text-align:right;"> 0.1437509 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age </td>
   <td style="text-align:left;"> 25-34 </td>
   <td style="text-align:right;"> 0.1515670 </td>
   <td style="text-align:right;"> 0.1515679 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age </td>
   <td style="text-align:left;"> 35-54 </td>
   <td style="text-align:right;"> 0.3299868 </td>
   <td style="text-align:right;"> 0.3299894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age </td>
   <td style="text-align:left;"> 55+ </td>
   <td style="text-align:right;"> 0.3246185 </td>
   <td style="text-align:right;"> 0.3246167 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sex </td>
   <td style="text-align:left;"> Male </td>
   <td style="text-align:right;"> 0.4890179 </td>
   <td style="text-align:right;"> 0.4890137 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sex </td>
   <td style="text-align:left;"> Female </td>
   <td style="text-align:right;"> 0.5109821 </td>
   <td style="text-align:right;"> 0.5109863 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Education </td>
   <td style="text-align:left;"> Less than secondary </td>
   <td style="text-align:right;"> 0.1853487 </td>
   <td style="text-align:right;"> 0.1853412 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Education </td>
   <td style="text-align:left;"> Secondary school grad </td>
   <td style="text-align:right;"> 0.2000918 </td>
   <td style="text-align:right;"> 0.2000976 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Education </td>
   <td style="text-align:left;"> Some post-secondary </td>
   <td style="text-align:right;"> 0.0515635 </td>
   <td style="text-align:right;"> 0.0515641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Education </td>
   <td style="text-align:left;"> Post-secondary grad </td>
   <td style="text-align:right;"> 0.5489343 </td>
   <td style="text-align:right;"> 0.5489352 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Education </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0140618 </td>
   <td style="text-align:right;"> 0.0140620 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> No income </td>
   <td style="text-align:right;"> 0.0250418 </td>
   <td style="text-align:right;"> 0.0250420 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> LESS THAN $5,000 </td>
   <td style="text-align:right;"> 0.0282772 </td>
   <td style="text-align:right;"> 0.0282781 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $5,000 TO $9,999 </td>
   <td style="text-align:right;"> 0.0424744 </td>
   <td style="text-align:right;"> 0.0424748 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $10,000 TO $14,999 </td>
   <td style="text-align:right;"> 0.0641120 </td>
   <td style="text-align:right;"> 0.0641134 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $15,000 TO $19,999 </td>
   <td style="text-align:right;"> 0.0512046 </td>
   <td style="text-align:right;"> 0.0512049 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $20,000 TO $29,999 </td>
   <td style="text-align:right;"> 0.1010808 </td>
   <td style="text-align:right;"> 0.1010808 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $30,000 TO $39,999 </td>
   <td style="text-align:right;"> 0.0898432 </td>
   <td style="text-align:right;"> 0.0898450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $40,000 TO $49,999 </td>
   <td style="text-align:right;"> 0.0876461 </td>
   <td style="text-align:right;"> 0.0876457 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $50,000 TO $59,999 </td>
   <td style="text-align:right;"> 0.0672547 </td>
   <td style="text-align:right;"> 0.0672537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $60,000 TO $79,999 </td>
   <td style="text-align:right;"> 0.0838845 </td>
   <td style="text-align:right;"> 0.0838832 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $80,000 TO 99,999 </td>
   <td style="text-align:right;"> 0.0535977 </td>
   <td style="text-align:right;"> 0.0535978 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> $100k or more </td>
   <td style="text-align:right;"> 0.0701623 </td>
   <td style="text-align:right;"> 0.0701628 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 0.0789694 </td>
   <td style="text-align:right;"> 0.0789660 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Personal income </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.1564515 </td>
   <td style="text-align:right;"> 0.1564520 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> No income </td>
   <td style="text-align:right;"> 0.0034083 </td>
   <td style="text-align:right;"> 0.0034083 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> LESS THAN $5,000 </td>
   <td style="text-align:right;"> 0.0039214 </td>
   <td style="text-align:right;"> 0.0039213 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $5,000 TO $9,999 </td>
   <td style="text-align:right;"> 0.0082869 </td>
   <td style="text-align:right;"> 0.0082870 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $10,000 TO $14,999 </td>
   <td style="text-align:right;"> 0.0226351 </td>
   <td style="text-align:right;"> 0.0226357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $15,000 TO $19,999 </td>
   <td style="text-align:right;"> 0.0321571 </td>
   <td style="text-align:right;"> 0.0321582 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $20,000 TO $29,999 </td>
   <td style="text-align:right;"> 0.0842706 </td>
   <td style="text-align:right;"> 0.0842693 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $30,000 TO $39,999 </td>
   <td style="text-align:right;"> 0.0850972 </td>
   <td style="text-align:right;"> 0.0850999 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $40,000 TO $49,999 </td>
   <td style="text-align:right;"> 0.0896067 </td>
   <td style="text-align:right;"> 0.0896063 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $50,000 TO $59,999 </td>
   <td style="text-align:right;"> 0.0810646 </td>
   <td style="text-align:right;"> 0.0810649 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $60,000 TO $79,999 </td>
   <td style="text-align:right;"> 0.1411437 </td>
   <td style="text-align:right;"> 0.1411425 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $80,000 TO 99,999 </td>
   <td style="text-align:right;"> 0.1129619 </td>
   <td style="text-align:right;"> 0.1129603 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> $100k or more </td>
   <td style="text-align:right;"> 0.3331578 </td>
   <td style="text-align:right;"> 0.3331580 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Household income </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0022887 </td>
   <td style="text-align:right;"> 0.0022885 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> White </td>
   <td style="text-align:right;"> 0.6898763 </td>
   <td style="text-align:right;"> 0.6898760 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Black </td>
   <td style="text-align:right;"> 0.0381555 </td>
   <td style="text-align:right;"> 0.0381564 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Korean </td>
   <td style="text-align:right;"> 0.0037696 </td>
   <td style="text-align:right;"> 0.0037698 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Filipino </td>
   <td style="text-align:right;"> 0.0238718 </td>
   <td style="text-align:right;"> 0.0238720 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Japanese </td>
   <td style="text-align:right;"> 0.0019801 </td>
   <td style="text-align:right;"> 0.0019799 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Chinese </td>
   <td style="text-align:right;"> 0.0425515 </td>
   <td style="text-align:right;"> 0.0425516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> South Asian </td>
   <td style="text-align:right;"> 0.0653306 </td>
   <td style="text-align:right;"> 0.0653301 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> South East Asian </td>
   <td style="text-align:right;"> 0.0121126 </td>
   <td style="text-align:right;"> 0.0121130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Arab </td>
   <td style="text-align:right;"> 0.0092517 </td>
   <td style="text-align:right;"> 0.0092518 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> West Asian </td>
   <td style="text-align:right;"> 0.0106736 </td>
   <td style="text-align:right;"> 0.0106735 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Latin American </td>
   <td style="text-align:right;"> 0.0174512 </td>
   <td style="text-align:right;"> 0.0174516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Other Racial or Cultural Origin </td>
   <td style="text-align:right;"> 0.0140982 </td>
   <td style="text-align:right;"> 0.0140985 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> Multiple Racial/Cultural Origins </td>
   <td style="text-align:right;"> 0.0164220 </td>
   <td style="text-align:right;"> 0.0164218 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 0.0241613 </td>
   <td style="text-align:right;"> 0.0241616 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ethnicity </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0302941 </td>
   <td style="text-align:right;"> 0.0302924 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Married </td>
   <td style="text-align:right;"> 0.5032409 </td>
   <td style="text-align:right;"> 0.5032429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Common-Law </td>
   <td style="text-align:right;"> 0.0660569 </td>
   <td style="text-align:right;"> 0.0660579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Widowed </td>
   <td style="text-align:right;"> 0.0458991 </td>
   <td style="text-align:right;"> 0.0458977 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Separated </td>
   <td style="text-align:right;"> 0.0258048 </td>
   <td style="text-align:right;"> 0.0258058 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Divorced </td>
   <td style="text-align:right;"> 0.0487190 </td>
   <td style="text-align:right;"> 0.0487166 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Single/Never Married </td>
   <td style="text-align:right;"> 0.3075746 </td>
   <td style="text-align:right;"> 0.3075741 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> DK </td>
   <td style="text-align:right;"> 0.0010358 </td>
   <td style="text-align:right;"> 0.0010359 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> Refusal </td>
   <td style="text-align:right;"> 0.0012782 </td>
   <td style="text-align:right;"> 0.0012783 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Marital status </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0003908 </td>
   <td style="text-align:right;"> 0.0003908 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Geography </td>
   <td style="text-align:left;"> Population Centre </td>
   <td style="text-align:right;"> 0.8429496 </td>
   <td style="text-align:right;"> 0.8429489 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Geography </td>
   <td style="text-align:left;"> Rural </td>
   <td style="text-align:right;"> 0.1570504 </td>
   <td style="text-align:right;"> 0.1570511 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> less than once a month </td>
   <td style="text-align:right;"> 0.1687812 </td>
   <td style="text-align:right;"> 0.1687812 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> once a month </td>
   <td style="text-align:right;"> 0.0745744 </td>
   <td style="text-align:right;"> 0.0745746 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> 2 to 3 times a month </td>
   <td style="text-align:right;"> 0.1006304 </td>
   <td style="text-align:right;"> 0.1006314 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> once a week </td>
   <td style="text-align:right;"> 0.1110777 </td>
   <td style="text-align:right;"> 0.1110780 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> 2 to 3 times a week </td>
   <td style="text-align:right;"> 0.1517592 </td>
   <td style="text-align:right;"> 0.1517594 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> 4 to 6 times a week </td>
   <td style="text-align:right;"> 0.0452207 </td>
   <td style="text-align:right;"> 0.0452212 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> every day </td>
   <td style="text-align:right;"> 0.0665190 </td>
   <td style="text-align:right;"> 0.0665203 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 0.2623282 </td>
   <td style="text-align:right;"> 0.2623253 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> DK </td>
   <td style="text-align:right;"> 0.0019352 </td>
   <td style="text-align:right;"> 0.0019353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> Refusal </td>
   <td style="text-align:right;"> 0.0006117 </td>
   <td style="text-align:right;"> 0.0006120 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alcohol drinking frequency </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0165623 </td>
   <td style="text-align:right;"> 0.0165612 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Type of alcohol drinker </td>
   <td style="text-align:left;"> Regular drinker </td>
   <td style="text-align:right;"> 0.5497813 </td>
   <td style="text-align:right;"> 0.5497851 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Type of alcohol drinker </td>
   <td style="text-align:left;"> Occasional drinker </td>
   <td style="text-align:right;"> 0.1687812 </td>
   <td style="text-align:right;"> 0.1687812 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Type of alcohol drinker </td>
   <td style="text-align:left;"> Did not drink in the last 12 months </td>
   <td style="text-align:right;"> 0.2623282 </td>
   <td style="text-align:right;"> 0.2623253 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Type of alcohol drinker </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0191092 </td>
   <td style="text-align:right;"> 0.0191084 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Heavy drinking </td>
   <td style="text-align:left;"> Heavy drinking in the past year at least once per month </td>
   <td style="text-align:right;"> 0.1628624 </td>
   <td style="text-align:right;"> 0.1628633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Heavy drinking </td>
   <td style="text-align:left;"> Did not heavy drink in the past year </td>
   <td style="text-align:right;"> 0.5536088 </td>
   <td style="text-align:right;"> 0.5536124 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Heavy drinking </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 0.2623282 </td>
   <td style="text-align:right;"> 0.2623253 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Heavy drinking </td>
   <td style="text-align:left;"> DK </td>
   <td style="text-align:right;"> 0.0027260 </td>
   <td style="text-align:right;"> 0.0027256 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Heavy drinking </td>
   <td style="text-align:left;"> Refusal </td>
   <td style="text-align:right;"> 0.0010243 </td>
   <td style="text-align:right;"> 0.0010245 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Heavy drinking </td>
   <td style="text-align:left;"> NS </td>
   <td style="text-align:right;"> 0.0174501 </td>
   <td style="text-align:right;"> 0.0174490 </td>
  </tr>
</tbody>
</table>

<img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-1.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-2.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-3.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-4.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-5.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-6.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-7.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-8.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-9.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-10.png" width="672" /><img src="cchs_data_trans_prob_sim_files/figure-html/comparison_barplots-11.png" width="672" />

# Transition probabilities

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Age </th>
   <th style="text-align:left;"> Sex </th>
   <th style="text-align:right;"> Alcoholic death rates (per 100K) </th>
   <th style="text-align:right;"> Non-alcoholic death rates (per 100K) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 16-24 </td>
   <td style="text-align:left;"> Male </td>
   <td style="text-align:right;"> 10.6 </td>
   <td style="text-align:right;"> 43.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16-24 </td>
   <td style="text-align:left;"> Female </td>
   <td style="text-align:right;"> 2.1 </td>
   <td style="text-align:right;"> 26.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25-34 </td>
   <td style="text-align:left;"> Male </td>
   <td style="text-align:right;"> 16.6 </td>
   <td style="text-align:right;"> 58.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25-34 </td>
   <td style="text-align:left;"> Female </td>
   <td style="text-align:right;"> 4.4 </td>
   <td style="text-align:right;"> 36.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35-54 </td>
   <td style="text-align:left;"> Male </td>
   <td style="text-align:right;"> 41.4 </td>
   <td style="text-align:right;"> 192.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35-54 </td>
   <td style="text-align:left;"> Female </td>
   <td style="text-align:right;"> 20.4 </td>
   <td style="text-align:right;"> 148.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55+ </td>
   <td style="text-align:left;"> Male </td>
   <td style="text-align:right;"> 72.1 </td>
   <td style="text-align:right;"> 2741.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55+ </td>
   <td style="text-align:left;"> Female </td>
   <td style="text-align:right;"> 28.9 </td>
   <td style="text-align:right;"> 3083.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total </td>
   <td style="text-align:left;"> Male </td>
   <td style="text-align:right;"> 42.8 </td>
   <td style="text-align:right;"> 1000.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total </td>
   <td style="text-align:left;"> Female </td>
   <td style="text-align:right;"> 17.8 </td>
   <td style="text-align:right;"> 1082.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> All </td>
   <td style="text-align:left;"> All </td>
   <td style="text-align:right;"> 30.3 </td>
   <td style="text-align:right;"> 1041.0 </td>
  </tr>
</tbody>
</table>




# Combine alcoholic and non-alcholic attributable transition probabilities

## Poor method


```{=html}
<div id="htmlwidget-9b5348d04508dd795bb3" style="width:672px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-9b5348d04508dd795bb3">{"x":{"diagram":"digraph flowchart {\n      # node definitions with substituted label text\n      graph [layout = dot, rankdir = LR]\n      node [fontname = Helvetica, shape = rectangle, style=filled]        \n      tab1 [label = \"Alcoholic attr. mort. rate\", color = orange]\n      tab2 [label = \"Non-alcoholic attr. mort. rate\", color = grey]\n      tab3 [label = \"Sample either \n  alcoholic or non-alcoholic \n  rate based on \n % mort. alcohol-attr.\", shape = diamond]\n      tab4 [label = \"Monte Carlo \n simulation on the \n selected mort. rate\", shape = diamond]\n      tab5 [label = \"Alcoholic attr. mort.\", color = pink]\n      \n      # edge definitions with the node IDs\n      tab1 -> tab3;\n      tab2 -> tab3;\n      tab3 -> tab4;\n      tab4 -> tab5\n      }\n      ","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>
```

This make not be correct, since it is possible that a person die from both alcoholic and non-alcoholic causes. Also not sure they alcoholic and non-alcoholic attributable rates are exclusive of each other.

# Better method


```{=html}
<div id="htmlwidget-4bea38b404cdd4164b32" style="width:672px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-4bea38b404cdd4164b32">{"x":{"diagram":"digraph flowchart {\n      # node definitions with substituted label text\n      graph [layout = dot, rankdir = LR]\n      node [fontname = Helvetica, shape = rectangle, style=filled]        \n      tab1 [label = \"Alcoholic att. mort. rate\", color = orange]\n      tab2 [label = \"Non-alcoholic att. mort. rate\", color = grey]\n      tab3 [label = \"Monte Carlo \n simulation\", shape = diamond, color = orange]\n      tab4 [label = \"Monte Carlo \n simulation\", shape = diamond, color = grey]\n      tab5 [label = \"Combine Mone Carlo\n results (summation)\", shape = diamond]\n      tab6 [label = \"Alcoholic att. mort.\", color = pink]\n      \n      # edge definitions with the node IDs\n      tab1 -> tab3;\n      tab2 -> tab4;\n      tab3 -> tab5;\n      tab4 -> tab5;\n      tab5 -> tab6\n      }\n\n      ","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>
```

This maybe a better solution. However, if alcoholic and non-alcoholic attributable rates are NOT exclusive of each other, then the mortality maybe inflated.

<!-- Run 3 state model -->



# Run two-state model


```r
# Load microsimulation functions
source('sim_functions_2state.R')

# Set memory limit
memory.limit(size=20000)

# Run microsimulations
comp.time <- list()

for (i in seq(cchsDataList)[1:8]) {

  # Data subset
  cchsData <- cchsDataList[[i]]
  transProbSub <- transProb[i,]
  
  # Model input
  n.i   <- sum(round(cchsData$WTS_S))   # number of simulated individuals, since it is weighted sum it needs to be rounded to join back to CCHS data
  n.t   <- 30                    # time horizon, 30 cycles
  v.n   <- c("H","D")            # the model states: Healthy (H), Dead (D)
  n.s   <- length(v.n)           # the number of health states
  v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state 
  d.e <- 0.03                    # equal discounting of costs and QALYs by 3%
  samplev.m <- 1                 # argument use for samplev() function 
  
  # Transition probabilities (per cycle)
  p.HD1   <- transProbSub$transition_4_rate_100k/100000      # probability to die when healthy (alcoholic attributable)
  p.HD2   <- transProbSub$transition_3_rate_100k/100000      # probability to die when healthy (non-alcoholic attributable)
  p.AA    <- transProbSub$percent_alcohol_attributable       # % mortalities alcohol-attributable
  
  # QALY inputs
  u.H     <- 1              # utility when healthy 
  v.x     <- runif(n.i, 0.95, 1.05) # vector capturing individuals' effect modifier at baseline
  
  # START SIMULATION
  cat('\n# Stratum', paste(transProbSub[c('age','sex')], collapse=':'),', n =',n.i, '\n')      # display the stratum
  
  p = Sys.time()
  sim <- MicroSimV2(v.M_1, n.i, n.t, v.n, samplev.m, X = v.x, d.e, seed = 100)  # run for treatment
  comp.time[[i]] = Sys.time() - p
  
  print(comp.time[[i]])
  
  # SAVE SIMULATION RESULTS TO FILE
  file_path <- file.path(dataDir,paste0('sim_',i,'_2state_v2.rds'))
  saveRDS(sim, file_path)
  
  rm(sim)
}
```

# Results from two-state model for males ages 16-24











