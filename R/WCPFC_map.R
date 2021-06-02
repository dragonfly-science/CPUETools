WCPFC_map_quant <- function(data, quant, lat = 'lat5', lon = 'lon5', by_year = F, trans = 'I', labq = quant, uncert = NULL){

  pdat <- data %>%
    {if (by_year) group_by(.,!!lon,!!lat,yy) else group_by(.,!!lon,!!lat)} %>%
    summarise(!!labq := mean(!!quant,na.rm=T))


  if(!is.null(uncert)) {
    require(multiscales)

    g1 <- ggplot(pdat)  +
      {if (class(pdat)!='sf') geom_tile(aes(x=lon5+2.5,y=lat5+2.5,fill=zip(!!quant,uncert)))} +
      {if (class(pdat)=='sf') geom_sf(aes(fill=zip(!!quant,uncert)),data=ol)}+
      bivariate_scale("fill",
                      pal_vsup(values = cividis(8),max_desat = 0.1),
                      limits = list(c(0, 1), c(0,1)),
                      breaks = list(waiver(), waiver()),
                      name = c("Relative abundance", "Relative CV"),
                      labels = list(waiver(), waiver()),
                      trans = list(trans,'identity'),
                      guide = "colourfan"
      )
  } else {

    g1 <- ggplot(pdat)  +
      {if (class(pdat)!='sf') geom_tile(aes(x=lon5+2.5,y=lat5+2.5,fill=!!quant))} +
      {if (class(pdat)=='sf') geom_sf(aes(fill=!!quant),data=ol)}+
      scale_fill_viridis_c(option='E', trans = trans)
  }


  g1+   geom_sf(data=world) +
    geom_sf(data=wcpfc, col='seagreen2',fill=NA) +
    coord_sf(xlim=c(110,225),ylim=c(-50,0)) +
    geom_hline(yintercept = 0, linetype=2) +
    theme_cowplot() +
    theme(legend.position = c(0.08,0.4)
    ) +
    labs(y='Latitude', x='Longitude')

  if (by_year) g1 + facet_wrap(~yy, ncol = 3)
  print(g1)
}
