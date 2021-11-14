NZ_map_quant <- function(data,
                            quant=NULL,
                            lat = 'start_latitude',
                            lon = 'start_longitude',
                            xlim = NULL,
                            ylim = NULL,
                            facet = NULL,
                            trans = 'identity',
                            labq = quant,
                            fun='I',
                            uncert = NULL,
                            limits=NULL,
                            fun_uncert='sum',
                            uncert_lab = uncert,
                            adj=0.05){

  require(rlang)

  lat <- sym(lat)
  lon <- sym(lon)
  labq <- sym(labq)
  if (!is.null(uncert)) uncert_lab<- sym(uncert_lab)

  # browser()
  if(fun != 'I'){
    pdat <- data %>%
    {if (!is.null(facet)) { facets <- sym(facet); group_by(.,!!lon,!!lat,!!facets)} else group_by(.,!!lon,!!lat)} %>%
    {if (!is.null(uncert)){
      summarise(., !!labq := do.call(fun,list(eval(parse_expr(quant)))),
                !!uncert_lab := do.call(fun_uncert,list(eval(parse_expr(uncert)))))
    } else {
      summarise(., !!labq := do.call(fun,list(eval(parse_expr(quant)))))
    }
    }} else {
      pdat <- data %>%
        {if (!is.null(uncert)){
          mutate(., !!labq := eval(parse_expr(quant)),
                    !!uncert_lab := eval(parse_expr(uncert)))
        } else {
          mutate(., !!labq := eval(parse_expr(quant)))
        }
    }}
  #browser()

  if(!is.null(uncert)) {
    require(multiscales)

    #browser()
    g1 <- ggplot(pdat)  +
      {if (class(pdat)!='sf') geom_tile(aes(x=!!lon,y=!!lat,fill=zip(!!labq,!!uncert_lab), col=zip(!!labq,!!uncert_lab)))} +
      {if (class(pdat)=='sf') geom_sf(aes(fill=zip(!!labq,!!uncert_lab), col=zip(!!labq,!!uncert_lab)))}+
      bivariate_scale("fill",
                       pal_vsup(values = viridis::cividis(32),
                                unc_levels=6,max_desat = 0.1),
                      limits = list({if(!is.null(limits)) limits else c(signif(min(pdat[[labq]],na.rm=T)*(1-adj),2), signif(max(pdat[[labq]],na.rm=T)*(1+adj),2))},
                                    c(signif(min(pdat[[uncert_lab]],na.rm=T),digits = 2),max(signif(max(pdat[[uncert_lab]],na.rm=T)*1.05,digits = 2),0.2))),
                      breaks = list(waiver(), c(signif(min(pdat[[uncert_lab]],na.rm=T),digits = 2),max(signif(pdat[[uncert_lab]],digits = 2),0.2,na.rm=T))),
                      name = c(as_label(labq), as_label(uncert_lab)),
                      labels = list(waiver(), function(x) c('Low', 'High')),
                      trans = list(trans,'identity'),
                      guide = "colourfan" )+
      bivariate_scale("color",
                      pal_vsup(values = viridis::cividis(32),
                               unc_levels=6,max_desat = 0.1),
                      limits = list({if(!is.null(limits)) limits else c(signif(min(pdat[[labq]],na.rm=T)*(1-adj),2), signif(max(pdat[[labq]],na.rm=T)*(1+adj),2))},
                                    c(signif(min(pdat[[uncert_lab]],na.rm=T),digits = 2),max(signif(max(pdat[[uncert_lab]],na.rm=T)*1.05,digits = 2),0.2))),
                      breaks = list(waiver(), c(signif(min(pdat[[uncert_lab]],na.rm=T),digits = 2),max(signif(pdat[[uncert_lab]],digits = 2),0.2,na.rm=T))),
                      name = c(as_label(labq), as_label(uncert_lab)),
                      labels = list(waiver(), function(x) c('Low', 'High')),
                      trans = list(trans,'identity'),
                      guide = "none" )
  } else {
    #browser()
    g1 <- ggplot(pdat)  +
      {if (class(pdat)!='sf') geom_tile(aes(x=!!lon,y=!!lat,fill=!!labq))} +
      {if (class(pdat)=='sf') geom_sf(aes(fill=!!labq,col=!!labq))}+
      scale_fill_viridis_c(option='E', trans = trans, limits = limits)+
      scale_colour_viridis_c(option='E', trans = trans, limits = limits)
  }


  g1 <- g1+   geom_sf(data=nz, fill='green4',alpha=0.2) +
    geom_hline(yintercept = 0, linetype=2) +
    theme_cowplot(font_size = 12) +
    theme(legend.position = c(0.8,0.6),
          legend.text = element_text(size=8)
    ) +
    labs(y='Latitude', x='Longitude')

  if (!is.null(facet)) {
    if (facet=='yy') g1 <- g1 + facet_wrap(~yy, ncol = floor(sqrt(length(unique(pdat$yy))))) + theme_cowplot(font_size = 12)
    if (facet=='flag_id') g1 <- g1 + facet_wrap(~flag_id, ncol = floor(sqrt(length(unique(pdat$flag_id))))) + theme_cowplot(font_size = 12)
    if (facet=='vessel_id') g1 <- g1 + facet_wrap(~vessel_id, ncol = floor(sqrt(length(unique(pdat$vessel_id))))) + theme_cowplot(font_size = 12)
    if (facet=='mm') g1 <- g1 + facet_wrap(~mm, ncol = floor(sqrt(length(unique(pdat$mm))))) + theme_cowplot(font_size = 12)
    if (facet=='Fleet') g1 <- g1 + facet_wrap(~Fleet, nrow = floor(sqrt(length(unique(pdat$Fleet))))) + theme_cowplot(font_size = 12)
  }
  g1 + coord_sf(xlim = xlim, ylim=ylim)
}
