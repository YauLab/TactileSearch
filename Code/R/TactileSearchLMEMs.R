
setwd('/Volumes/bcm-neuro-yau/data/Liz/Code/R/PBS/')
source('/Volumes/bcm-neuro-yau/data/Liz/Code/R/PBS/basic_plotters.R')

require(lme4)
require(lmerTest)
require(readxl)

# Sync data ---------------------------------------------------------------
dsSync <- readxl::read_excel('SearchData_Sync.xls')
# Rename so order of fixed effects is desired:
dsSync$TPresent <- chartr('NO','bAb',dsSync$TPresent)
dsSync$TPresent <- chartr('YES','aPr',dsSync$TPresent)
dsSync$DSync <- chartr('NO','bAb',dsSync$DSync)
dsSync$DSync <- chartr('YES','aPr',dsSync$DSync)
# Turn response times into ms:
dsSync$RT <- dsSync$RT * 1000
# LME:
modSync.lmm <- lmer(RT ~ SetSize*TPresent*TFreq + (1|Subject), data=dsSync, REML=FALSE)
summary(modSync.lmm)

# PhSc data ---------------------------------------------------------------
dsPhSc <- readxl::read_excel('SearchData_PhSc.xls')
# Rename so order of fixed effects is desired:
dsPhSc$TPresent <- chartr('NO','bAb',dsPhSc$TPresent)
dsPhSc$TPresent <- chartr('YES','aPr',dsPhSc$TPresent)
dsPhSc$DSync <- chartr('NO','bAb',dsPhSc$DSync)
dsPhSc$DSync <- chartr('YES','aPr',dsPhSc$DSync)
# Turn response times into ms:
dsPhSc$RT <- dsPhSc$RT * 1000
# LME:
modPhSc.lmm <- lmer(RT ~ SetSize*TPresent*TFreq + (1|Subject), data=dsPhSc, REML=FALSE)
summary(modPhSc.lmm)

# Combined data -----------------------------------------------------------
dsComb <- readxl::read_excel('SearchDataNew.xls')
# Rename so order of fixed effects is desired:
dsComb$TPresent <- chartr('NO','bAb',dsComb$TPresent)
dsComb$TPresent <- chartr('YES','aPr',dsComb$TPresent)
dsComb$DSync <- chartr('NO','bAb',dsComb$DSync)
dsComb$DSync <- chartr('YES','aPr',dsComb$DSync)
# Turn response times into ms:
dsComb$RT <- dsComb$RT * 1000
# LME:
modComb.lmm <- lmer(RT ~ SetSize*TPresent*TFreq*DSync + (1|Subject), data=dsComb, REML=FALSE)
summary(modComb.lmm)
# Predicted data:
dsComb$yhat.lmm <- predict(modComb.lmm)



