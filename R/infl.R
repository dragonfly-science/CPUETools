
plot_CEs <- function(model,
                     fx,
                     idx = 'yy',
                     resp = NULL,
                     varlab = NULL,
                     trans_eff=NULL,
                     trans_resp=NULL,
                     limits=NULL,
                     dpar=NULL){

  if(length(varlab)==2) {

    fx1<- strsplit(fx, ':')[[1]][1]
    fx2<- strsplit(fx, ':')[[1]][2]

    if(!is.null(attr(model$data[[fx2]],'scaled:scale'))) {
      cdat <- model$data[[fx2]]*attr(model$data[[fx2]],'scaled:scale')+attr(model$data[[fx2]],'scaled:center')
      cond <- c(quantile(cdat,0), mean(cdat), quantile(cdat,1))
      cond <- (cond-attr(model$data[[fx2]],'scaled:center'))/attr(model$data[[fx2]],'scaled:scale')
      conds <- data.frame(cond)
      colnames(conds) <- fx2
      CE <- conditional_effects(model, effects = x, int_conditions = conds, plot=F,re_formula=NA, dpar = dpar)
    } else {
      CE <- conditional_effects(model, effects = fx, plot=F,re_formula=NA, dpar = dpar)
    }
  }  else {
    CE <- conditional_effects(model, effects = fx, plot=F,re_formula=NA, dpar = dpar)
  }

  if(length(varlab)==2) {

    if(!is.null(attr(model$data[[fx1]],'scaled:center'))) {
      CE[[fx]]['effect1__'] <- CE[[fx]]['effect1__']*attr(model$data[[fx1]],'scaled:scale')+attr(model$data[[fx1]],'scaled:center')
      CE[[fx]][fx1] <- CE[[fx]][fx1]*attr(model$data[[fx1]],'scaled:scale')+attr(model$data[[fx1]],'scaled:center')
    }
    if(!is.null(attr(model$data[[fx2]],'scaled:center'))) {
      CE[[fx]]['effect2__'] <-    factor(signif(unlist(CE[[fx]][fx2]*attr(model$data[[fx2]],'scaled:scale')+attr(model$data[[fx2]],'scaled:center')),2))
      CE[[fx]][fx2] <- CE[[fx]][fx2]*attr(model$data[[fx2]],'scaled:scale')+attr(model$data[[fx2]],'scaled:center')

    }
  } else {
    if(!is.null(attr(model$data[[fx]],'scaled:center'))) {

      CE[[fx]][fx] <- CE[[fx]][fx]*attr(model$data[[fx]],'scaled:scale')+attr(model$data[[fx]],'scaled:center')
      CE[[fx]]['effect1__'] <- CE[[fx]]['effect1__']*attr(model$data[[fx]],'scaled:scale')+attr(model$data[[fx]],'scaled:center')
    }
  }

  if(!is.null(resp)) attr(CE[[fx]], "response") <- resp
  if(!is.null(varlab)) attr(CE[[fx]], "effects") <- varlab

  if(!is.null(trans_eff)) CE[[fx]]['effect1__'] <- do.call(trans_eff, list(CE[[fx]]['effect1__']))
  #browser()
  if(!is.null(trans_resp)) CE[[fx]][c('estimate__','lower__','upper__')] <- do.call(trans_resp, list(CE[[fx]][c('estimate__','lower__','upper__')]))
  if(!is.null(limits)) CE[[fx]] %<>%  filter(effect1__>=limits[1], effect1__<=limits[2])

  CE[[fx]]$upper__[CE[[fx]]$upper__ > 5*max(CE[[fx]]$estimate__) ] <- 5*max(CE[[fx]]$estimate__)

  if(!is.null(CE[[fx]]$effect2__)) {

    g <- ggplot(CE[[fx]], aes(x=effect1__, y= estimate__)) +
      xlab(attr(CE[[fx]],'effects')[1]) +
      ylab(attr(CE[[fx]],'response')[1]) +
      scale_color_viridis_d(attr(CE[[fx]],'effects')[2], option='E')

    if(is.factor(CE[[fx]]$effect1__)){
      g <- g + geom_pointrange(aes(ymin = lower__, ymax=upper__, col=effect2__))
    } else {
      g <- g + geom_line(aes(col=effect2__))+
        geom_ribbon(aes(ymin = lower__, ymax=upper__, col=effect2__), alpha=0.2)
    }

  } else {

    g <- ggplot(CE[[fx]],aes(x=effect1__, y= estimate__)) +
      xlab(attr(CE[[fx]],'effects')[1]) +
      ylab(attr(CE[[fx]],'response')[1])

    if(is.factor(CE[[fx]]$effect1__)){
      g <- g + geom_pointrange(aes(ymin = lower__, ymax=upper__))
    } else {
      g <- g + geom_line()+
        geom_ribbon(aes(ymin = lower__, ymax=upper__), alpha=0.2)
    }

  }
  return(g)
}


get_rand_Eff <- function(pred_data, fx, bmod){

  cpuedf <- pred_data %>%
    mutate(iid=paste0('r_',fx,'[', !!sym(fx), ',Intercept]'))
  ps <- colMeans(as_draws_df(bmod, variable =paste0('r_',fx,'\\['),regex = T) %>% select(-starts_with('.')))
  if(any(grepl('hu',bmod$formula))) pss <- try(colMeans(as_draws_df(bmod, variable =paste0('r_',fx,'__hu\\['),regex = T) %>% select(-starts_with('.')))) else pss <- try(a+b, silent=T)
  if(!class(pss) == 'try-error') {
    pn <- names(ps)
    #browser()
    ps <- log((1-inv_logit(pss))*exp(colMeans(as_draws_df(bmod,variable = "b_Intercept") %>% select(-starts_with('.')))+ps))
    ps = ps - mean(ps)
    names(ps) <- pn
  }

  cpuedf$Eff <- ps[gsub(' ','.',cpuedf$iid)]
  cpuedf
}

plot_infl <- function(fx,
                      idx = 'yy',
                      lab,
                      bmod,
                      plot_ylabs=F,
                      resp_lab = "CPUE (n/100 hooks)",
                      pred_data = NULL,
                      inord = F,
                      trans=NULL,
                      xlabs='Year'){

  if(is.null(pred_data)) pred_data <- bmod$data
  theme_set(theme_cowplot(font_size=12))
  blues <- colorRampPalette(c('dodgerblue4', 'royalblue', 'lightsteelblue'))
  reds <- colorRampPalette(c('peachpuff', 'tomato', 'tomato4'))
  influpal <- c(blues(20),reds(20))
  wcovar <- 'cluster'
  wvarint <- c(cluster='b_cluster2')
  relpal <- c('royalblue', 'orangered4')


  num <- is.numeric(pred_data[[fx]])
  rand <- fx %in% names(ranef(bmod))

  if(rand){

    cpuedf <- get_rand_Eff(pred_data, fx, bmod)

    if(!is.null(pred_data))  {
      Effmean <- get_rand_Eff(bmod$data, fx, bmod) %>% pull(Eff) %>% mean(na.rm=T)
    } else {
      Effmean = cpuedf %>% pull(Eff) %>% mean(na.rm=T)
    }

  } else if (!num){
    cpuedf <- pred_data

    ps <- colMeans(posterior_samples(bmod, pars=c('b_Intercept',paste0('b_',fx))))
    ps[2:length(ps)] <- ps[1] + ps[2:length(ps)]
    cpuedf$Eff <- ps[as.numeric(as.factor(cpuedf[[fx]]))]
    Effmean = cpuedf %>% pull(Eff) %>% mean(na.rm=T)
  } else {
    #browser()
    ce <- conditional_effects(bmod,fx,plot=F,resolution = 10,re_formula = NA)
    fxc <- unique(ce[[fx]]$effect1__)
    udiff <- diff(fxc)[1]/2
    ll <- last(fxc)

    cecond <- ce[[fx]] %>%
      mutate(conds = as.numeric(as.factor(effect1__))) %>%
      select(conds, Eff = estimate__)

#browser()
    labs = c(fxc - udiff, ll + udiff)
    if(!is.null(attr(pred_data[[fx]],'scaled:scale')))  labs <- labs*attr(pred_data[[fx]],'scaled:scale')+attr(pred_data[[fx]],'scaled:center')
    if(!is.null(trans)) labs <- do.call(trans, list(labs))
    labs <- paste('<=',round(labs,2))

    cpuedf <- pred_data %>%
      mutate(cond_fx = cut(!!sym(fx), breaks=c(-Inf,fxc - udiff, ll + udiff), labels=labs),
             conds = as.numeric(cond_fx)) %>%
      inner_join(cecond) %>%
      mutate(Eff = log(Eff/gmean(Eff)))

    if(!is.null(pred_data))  {
      Effmean = bmod$data %>%
        mutate(cond_fx = cut(!!sym(fx), breaks=c(fxc - udiff, ll + udiff)),
               conds = as.numeric(cond_fx)) %>%
        inner_join(cecond) %>%
        mutate(Eff = log(Eff/gmean(Eff))) %>% pull(Eff) %>% mean(na.rm=T)
    } else {
      Effmean = cpuedf %>% pull(Eff) %>% mean(na.rm=T)
    }

    cpuedf %<>% mutate(!!sym(fx) := cond_fx)

  }
  #browser()
  infldf  <- cpuedf %>%
    mutate(centr=Eff-Effmean) %>%
    group_by(!!sym(idx)) %>%
    summarize(neff=n(),
              suminf=sum(centr,na.rm=T),
              avinf=suminf/neff,
              Influ=exp(avinf)) %>%
    mutate(rel=ifelse(Influ<1, 'Negative', 'Positive'))

  g1 <- ggplot(infldf, aes(!!sym(idx), Influ,  colour=rel, group=NA)) +
    geom_hline(yintercept=1, linetype='dotted') +
    geom_line(colour='grey50') +
    geom_point(size=4) + ylab('Influence\n') + xlab('') +
    #    theme(legend.position='none') +
    scale_colour_manual('Influence',values=rev(relpal)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = ggplot2::margin(t = 10, r=0.5, b = 0, l = 1))



  covardf <- cpuedf %>% group_by(!!sym(idx), !!sym(fx)) %>%
    summarize(num=n(),
              Eff = exp(unique(Eff)),
              .groups = "drop")

  pclims <- max(abs(covardf$Eff))
  if(!is.null(attr(covardf[[fx]],'scaled:center'))) covardf[fx] <- covardf[[fx]]*attr(covardf[[fx]],'scaled:scale')+attr(covardf[[fx]],'scaled:center')

  if((!num | rand) & !inord){
    covardf %<>% ungroup() %>% arrange(Eff)
  } else {
    covardf %<>% ungroup() %>% arrange(!!sym(fx))
  }
  covardf %<>% mutate(fx = factor(!!sym(fx), levels=unique(!!sym(fx)), ordered=TRUE))


  g2 <- ggplot(covardf, aes(!!sym(idx), fx, size=num, colour=Eff, fill=Eff)) +
    geom_point(shape=21, colour='grey50') +
    scale_size('Events',range=c(1,10)) +
    ylab(lab) + xlab(xlabs) +
    theme(axis.text.x = element_text(angle=45,hjust=1))

  if(!rand) g2 <- g2 +
    scale_colour_gradientn(resp_lab, colours=rev(influpal)) +
    scale_fill_gradientn(resp_lab,colours=rev(influpal))


  if(rand | num) g2 <- g2 +
    scale_colour_gradient2('Multiplier', low = 'tomato4', high = 'dodgerblue4',midpoint = 1) +
    scale_fill_gradient2('Multiplier', low = 'tomato4', high = 'dodgerblue4',midpoint = 1)

  if(!plot_ylabs) g2 = g2+ theme(axis.text.y = element_blank())

  require(patchwork)
  g1/g2 + patchwork::plot_layout(guides = 'collect')
}
