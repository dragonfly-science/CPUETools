WCPFC_map_quant <- function(data,
                            quant=NULL,
                            lat = 'lat5',
                            lon = 'lon5',
                            facet = NULL,
                            trans = 'identity',
                            labq = quant,
                            fun='sum',
                            uncert = NULL,
                            limits=NULL,
                            fun_uncert='sum',
                            uncert_lab = uncert,
                            adj=0.05,
                            y1=-50,
                            y2=50){

  require(rlang)

  lat <- sym(lat)
  lon <- sym(lon)
  labq <- sym(labq)
  if (!is.null(uncert)) uncert_lab<- sym(uncert_lab)

  # browser()
  pdat <- data %>%
    {if (!is.null(facet)) { facets <- sym(facet); group_by(.,!!lon,!!lat,!!facets)} else group_by(.,!!lon,!!lat)} %>%
    {if (!is.null(uncert)){
      summarise(., !!labq := do.call(fun,list(eval(parse_expr(quant)))),
                !!uncert_lab := do.call(fun_uncert,list(eval(parse_expr(uncert)))))
    } else {
      summarise(., !!labq := do.call(fun,list(eval(parse_expr(quant)))))
    }
    }
  #browser()

  if(!is.null(uncert)) {
    require(multiscales)

    #browser()
    g1 <- ggplot(pdat)  +
      {if (any(class(pdat)!='sf')) geom_tile(aes(x=!!lon,y=!!lat,fill=zip(!!labq,!!uncert_lab)))} +
      {if (any(class(pdat)=='sf')) geom_sf(aes(fill=zip(!!labq,!!uncert_lab)))}+
      bivariate_scale("fill",
                      pal_vsup(values = viridis::viridis(32),
                               unc_levels=6,max_desat = 0.1),
                      limits = list(c(signif(min(pdat[labq],na.rm=T)*(1-adj),2), signif(max(pdat[labq],na.rm=T)*(1+adj),2)),
                                    c(0,max(signif(max(pdat[uncert_lab],na.rm=T)*1.05,digits = 2),0.2))),
                      breaks = list(waiver(), c(signif(min(pdat[uncert_lab],na.rm=T),digits = 2)*0.02,max(signif(pdat[uncert_lab],digits = 2),0.2,na.rm=T)*0.98)),
                      name = c(as_label(labq), as_label(uncert_lab)),
                      labels = list(waiver(), function(x) c('Low', 'High')),
                      trans = list(trans,'identity'),
                      guide = "colourfan" )
  } else {
    #browser()
    g1 <- ggplot(pdat)  +
      {if (any(class(pdat)!='sf')) geom_tile(aes(x=!!lon,y=!!lat,fill=!!labq))} +
      {if (any(class(pdat)=='sf')) geom_sf(aes(fill=!!labq))}+
      scale_fill_viridis_c(option='E', trans = trans, limits = limits) +
      scale_colour_viridis_c(option='E', trans = trans, limits = limits)
  }


  g1 <- g1+   geom_sf(data=world) +
    geom_sf(data=wcpfc, col='seagreen2',fill=NA) +
    coord_sf(xlim=c(110,225),ylim=c(y1,y2)) +
    geom_hline(yintercept = 0, linetype=2) +
    theme_cowplot(font_size = 12) +
    theme(legend.position = c(0.06,0.2),
          legend.text = element_text(size=8)
    ) +
    labs(y='Latitude', x='Longitude')

  if (!is.null(facet)) {
    if (facet=='yy') g1 <- g1 + facet_wrap(~yy, ncol = floor(sqrt(length(unique(pdat$yy))))) + theme_cowplot(font_size = 12)
    if (facet=='flag_id') g1 <- g1 + facet_wrap(~flag_id, ncol = floor(sqrt(length(unique(pdat$flag_id))))) + theme_cowplot(font_size = 12)
    if (facet=='program_code') g1 <- g1 + facet_wrap(~program_code, ncol = floor(sqrt(length(unique(pdat$program_code))))) + theme_cowplot(font_size = 12)
    if (facet=='vessel_id') g1 <- g1 + facet_wrap(~vessel_id, ncol = floor(sqrt(length(unique(pdat$vessel_id))))) + theme_cowplot(font_size = 12)
    if (facet=='mm') g1 <- g1 + facet_wrap(~mm, ncol = floor(sqrt(length(unique(pdat$mm))))) + theme_cowplot(font_size = 12)
    if (facet=='Fleet') g1 <- g1 + facet_wrap(~Fleet, nrow = floor(sqrt(length(unique(pdat$Fleet))))) + theme_cowplot(font_size = 12)
  }
  g1
}
