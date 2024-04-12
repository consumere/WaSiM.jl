# s. Julia/R code:
# make markers
using RCall

R"
require(terra);require(magrittr)
lk='C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/wasim_pest/app/input'
setwd(lk)
ezg=rast('sinn100.ezg')
v0<-raster::values(ezg) %>% table() %>% names %>% noquote()
"

R"list(rbind(
  v0 %>% paste0('$set $dr', ., ' = ', '~ dr',.,'~'),
  v0 %>% paste0('$set $ki', ., ' = ', '~ ki',.,'~'),
  v0 %>% paste0('$set $kd', ., ' = ', '~ kd',.,'~'), 
  v0 %>% paste0('$set $kb', ., ' = ', '~ kb',.,'~'),
  v0 %>% paste0('$set $q0', ., ' = ', '~ q0',.,'~'),
  v0 %>% paste0('$set $sd', ., ' = ', '~ sd',.,'~')
)) %>%  unlist(., use.names = F) %>% clipr::write_clip()  
"
# $set $T0R = ~ t0r~	#T0R	temperature	limit	for	rain	(Grad	Celsius)
# $set $T0  = ~ t0~	#T0	temperature	limit	snow	melt

op()

function cb(x)
    return clipboard(x)
end

function pwc()
    
end
cb(wslpath())

zp(pwc)