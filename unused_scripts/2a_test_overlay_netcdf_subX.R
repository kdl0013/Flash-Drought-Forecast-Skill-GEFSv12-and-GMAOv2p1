#!/usr/bin/env R
rm(list=ls())

library(raster)




library(dplyr)
library(vroom)

library(ggmap)
library(mapdata)
library(ggplot2)
library(stringr)
library(ggpubr)
library(repr)
library(cowplot)
library(RColorBrewer)
library(grid)
require(viridis)
require(gridExtra)
require(corrplot)
require(Kendall)
require(dplyr)
require(zoo)


cpc_path <- ('~/Insync/OneDrive/NRT_CPC_Internship')
mask_path <- paste0(cpc_path,'/Data/CONUS_mask')
data_path <- paste0(cpc_path,'/Data/GMAO')


conus_mask <- raster(paste0(mask_path,'/NCA-LDAS_masks_SubX.nc4'),varname=USDM-HP_mask)
  


colfunc_wyr <<- colorRampPalette(c('white','yellow', 'red'))
colfunc_wyor <<- colorRampPalette(c('white','yellow', 'orange', 'red'))
colfunc_wlsofd <<- colorRampPalette(c('white','goldenrod3',
                                        'firebrick1','darkred'))
colfunc_wyfr <- colorRampPalette(c('white','yellow','firebrick1','darkred'))


#
#
######################################GRAPHS BEGIN
# Percent of Weeks in Flash Drought - Spatial Distribution -COMPLETED -------------------------------
#Other options explored 'eddi.D1', 'eddi.D3', 'smvi', 'smeddi'
# x11()

index_='eddi.D2'
spatial_3_def <- function(index_){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  directory_create <- ('Spatial_distribution')
  dir.create(directory_create)
  
  #need 1 EDDI file to overlay new cluster information on top of
  outpath = paste(box_, directory_create, sep = '/')
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #Open mask
  mask_file <- mask_wrapper(index_ = index_)
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)

  spatial_file = paste0(index_,'_weekly_FD_percentage_reviewer_suggestion.csv')
  
  fh <- read.csv(spatial_file, row.names =1, header =T)
  
  #Overlay cluster information on top of ascii file to not create a new map
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0){
      fh_ascii@data[row,1] = NA
    }else{
      fh_ascii@data[row,1] = fh[row,1]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_ = index_, fh_ascii = fh_ascii)
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  #Map of the CONUS
  baseData <- map_data('state')
  #Set title of ggplot
  title_name <- paste(toupper(index_),'\n Weekly FD Percentage Count', sep ='')
  #Set color pallette
  colfunc <<- colorRampPalette(c('yellow', 'red'))
  #Set breaks
  breaks_spatial <<- ceiling(seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), by =10))
  
  #Minimum week count is in D1, maximum week count is in D3
  if (index_ == 'eddi.D2'){
    d2_brek <<-  ceiling(seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), by =2))
  } else if (index_ == 'smpd'){
    smpd_brek <<-  ceiling(seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), by =2))
  }
  
  # #Test for discrete color scale
  # spatial_discrete <- 
  #   
  #   ggplot()+
  #   geom_tile(data = fh_raster,aes(x=long, y=lat, fill = factor(Obs)))+
  #   theme(legend.title = element_blank(), legend.position = 'none')+
  #   geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
  #   theme_bw()+
  #   coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
  #   labs(title=title_name, x="Longitude", y="Latitude", fill = "Weeks in\n Flash Drought")+ 
  #   xlab(label = "Longitude")+
  #   theme(plot.title = element_text(hjust = 0.5, size=15)) +
  #   scale_fill_brewer(na.value='white', breaks=c(0,30,60,90,120,150,180,210,240,270,300),
  #                     palette = 'YlOrRd', )
  # 
  # 
  # 
  # 
  # +
  #   scale_fill_gradientn(colors = colfunc(100),
  #                        breaks = breaks_spatial,
  #                        labels= as.character(breaks_spatial), 
  #                        na.value = 'white')+
  #   theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))
  # 
  # 

  spatial_final <-  ggplot()+
    geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
    theme(legend.title = element_blank(), legend.position = 'none')+
    geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    theme_bw()+
    coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    labs(title=title_name, x="Longitude", y="Latitude", fill = "% of Weeks in\n Flash Drought")+ 
    xlab(label = "Longitude")+
    theme(plot.title = element_text(hjust = 0.5, size=15)) +
    scale_fill_gradientn(colors = colfunc_wyfr(40),
                         breaks = breaks_spatial,
                         labels= as.character(breaks_spatial), 
                         na.value = 'white')+
    theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))
  
  spatial_final
  ggsave(filename = paste0(outpath,'/',toupper(index_),'_weekly_FD_percentage', '.tiff') , width = 8, height = 5,dpi = 300)
  
  return(spatial_final)
  
}


'Plot spatial distribution for all indexes,running the function
will save the graphs as well'
#Plot figures
#SMPD spatial
smpd <- spatial_3_def('smpd')
smpd_brek #breaks

brek_val = 2


smpd_fin1 <- smpd+ xlab(label = NULL) + labs(title=NULL) +ylab(label='SMPD')+
  scale_fill_gradientn(colors = colfunc_wyfr(300),
                       breaks = seq(min(smpd_brek),max(smpd_brek),by=2),
                       labels= as.character(seq(min(smpd_brek),max(smpd_brek),by=2)),
                       limits =c(min(smpd_brek),max(smpd_brek)),
                       na.value = 'white')+
  scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
  scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
  theme(plot.title = element_text(size=22),
        axis.title.y =element_text(size=22))

smpd_no_legend <- smpd+ xlab(label = NULL) + labs(title=NULL) +ylab(label='SMPD')+
  scale_fill_gradientn(colors = colfunc_wyfr(300),
                       breaks = NULL,
                       labels= NULL,
                       limits =c(min(smpd_brek),max(smpd_brek)),
                       na.value = 'white')+
  scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
  scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
  theme(plot.title = element_text(size=22),
        axis.title.y =element_text(size=22))



#EDDI spatial
edd2 <- spatial_3_def('eddi.D2')
d2_brek #breaks
d2_brek[1]=3
d2_brek[6]=15



#change to same scale as smpd brek_
# d2_brek <- smpd_brek

#Make eddi plot
edd_fin <- edd2 +  labs(title=NULL) + xlab(label = NULL) +ylab(label='EDFD')+
  scale_fill_gradientn(colors = colfunc_wyfr(1200),
                       breaks = seq(min(d2_brek),max(d2_brek),by=2),
                       labels= as.character(seq(min(d2_brek),max(d2_brek),by=2)),
                       limits =c(min(d2_brek),max(d2_brek)),
                       na.value = 'white')+
  scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
  scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
  theme(plot.title = element_text(size=22),
        axis.title.y =element_text(size=22)) +
  theme(legend.key.height = unit(1.4, 'in'), legend.key.width = unit(0.7, 'cm'), legend.text = element_text(size=12))

edd_fin
#Make eddi plot smaller legend
edd_small_leg <- edd2 +  labs(title=NULL) + xlab(label = NULL) +ylab(label='EDFD')+
  scale_fill_gradientn(colors = colfunc_wyfr(1200),
                       breaks = seq(min(d2_brek),max(d2_brek),by=2),
                       labels= as.character(seq(min(d2_brek),max(d2_brek),by=2)),
                       limits =c(min(d2_brek),max(d2_brek)),
                       na.value = 'white')+
  scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
  scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
  theme(plot.title = element_text(size=22),
        axis.title.y =element_text(size=22)) 


#grab legend, remove legend again after this
edfd_legend <- g_legend(edd_fin)

#remove legend from edfd
#Make eddi plot
edd_no_legend <- edd2 +  labs(title=NULL) + xlab(label = NULL) +ylab(label='EDFD')+
  scale_fill_gradientn(colors = colfunc_wyfr(300),
                       breaks = NULL,
                       labels= NULL,
                       limits =c(min(d2_brek),max(d2_brek)),
                       na.value = 'white')+
  scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
  scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
  theme(plot.title = element_text(size=22),
        axis.title.y =element_text(size=22)) +
  theme(legend.key.height = unit(4.9, 'in'), legend.key.width = unit(0.7, 'cm'))


# ggarrange(edd_fin,smpd_fin1, common.legend = T, legend ='right', ncol = 1, nrow=2)

# ggarrange(edd_fin,smpd_fin,common.legend = T,nrow=2,legend = c('right'))
ggdraw() +
  draw_plot(edd_no_legend, x = 0.0, y = .50, width = .83, height = .45) +
  draw_plot(smpd_no_legend, x = 0.0, y = 0, width = .83, height = .45)+
  draw_plot(edfd_legend, x=0.88, y= 0.15, width= .07,height=.7)+
  draw_text("A", x = 0.05, y = 0.97,face='bold',size=20)+
  draw_text("B", x = 0.05, y = 0.47,face='bold',size=20)


ggsave(filename = paste0(new_dir,'/Weekly_FD_percentage_edfd_smpd.tiff'),dpi = 300,
       height = 8.58, width = 8.17, units = 'in')



#Seperate the plots and save
ggarrange(edd_small_leg,smpd_fin1,common.legend = F, legend = 'right',
          ncol=1)

ggsave(filename = paste0(new_dir,'/Weekly_FD_percentage_edfd_smpd_split_legends.tiff'),dpi = 300,
       height = 10.1, width = 8.58, units = 'in')




# Spatial maps, case study 2012-COMPLETED -------------------------------------------


edfd_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/edfd_month_average.csv')
smpd_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/smpd_month_average.csv')

#test
# month=6
# year=2012
# file=edfd_file
# month_name='June'

spatial_map <- function(month, month_name, year, file, index_){
  #find index of correct date
  index_num <- which(file$month==month & file$year==year)
  subset <- file[index_num,3:ncol(file)]
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  directory_create <- ('Spatial_distribution')
  dir.create(directory_create)
  
  #need 1 EDDI file to overlay new cluster information on top of
  outpath = paste(box_, directory_create, sep = '/')
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #Open mask
  mask_file <- mask_wrapper(index_ = index_)
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  

  #Overlay cluster information on top of ascii file to not create a new map
  for (row in 1:ncol(subset)){
    if (mask_fh[row,1] == 0 || subset[row]==0){
      fh_ascii@data[row,1] = NA
      
    }else{
      fh_ascii@data[row,1] = subset[row]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_ = index_, fh_ascii = fh_ascii)
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  #Map of the CONUS
  baseData <- map_data('state')

  #Set color pallette
  colfunc_map <- colorRampPalette(c('yellow','red'))
  #Set breaks
  breaks_spatial <<- seq(min(0), max(fh_raster$Obs, na.rm = T), by =0.25)
  

  
  # #Test for discrete color scale
  # spatial_discrete <- 
  #   
  #   ggplot()+
  #   geom_tile(data = fh_raster,aes(x=long, y=lat, fill = factor(Obs)))+
  #   theme(legend.title = element_blank(), legend.position = 'none')+
  #   geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
  #   theme_bw()+
  #   coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
  #   labs(title=title_name, x="Longitude", y="Latitude", fill = "Weeks in\n Flash Drought")+ 
  #   xlab(label = "Longitude")+
  #   theme(plot.title = element_text(hjust = 0.5, size=15)) +
  #   scale_fill_brewer(na.value='white', breaks=c(0,30,60,90,120,150,180,210,240,270,300),
  #                     palette = 'YlOrRd', )
  # 
  # 
  # 
  # 
  # +
  #   scale_fill_gradientn(colors = colfunc(100),
  #                        breaks = breaks_spatial,
  #                        labels= as.character(breaks_spatial), 
  #                        na.value = 'white')+
  #   theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))
  # 
  # 
  

  spatial_final <-
    ggplot()+
    geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
    theme(legend.title = element_blank(), legend.position = 'none')+
    geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    theme_bw()+
    coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    labs(title=month_name, x="Longitude", y="Latitude", fill = "% of Time\n in Drought")+ 
    xlab(label = "Longitude")+
    theme(plot.title = element_text(hjust = 0.5, size=15)) +
    scale_fill_gradientn(colors = colfunc_map(100),
                         breaks = breaks_spatial,
                         labels= as.character(breaks_spatial), 
                         na.value = 'white')+
    theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))+
      scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
      scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
      theme(plot.title = element_text(size=16,face='bold'),
            axis.title.y =element_text(size=11))
    
    
    if (index_ == 'eddi.D2'){
      name_o <- 'EDFD'
    }else{
      name_o <- 'SMPD'
    }
    
    grob <- grobTree(textGrob(name_o, x=0.9,  y=0.05, hjust=0,
                              gp=gpar(col="black", fontsize=10)))
  
  outmap <- spatial_final+annotation_custom(grob)

  return(outmap)
 }

# mar_edfd <- spatial_map(3,'March',2012,edfd_file,'eddi.D2')
# mar_smpd <- spatial_map(3,'March',2012,smpd_file,'smpd')

#Only do this analysis for reviewers comments
apr_edfd <- spatial_map(4,'April',2012,edfd_file,'eddi.D2')
apr_smpd <- spatial_map(4,'April',2012,smpd_file,'smpd')

may_edfd <- spatial_map(5,'May',2012,edfd_file,'eddi.D2')
may_smpd <- spatial_map(5,'May',2012,smpd_file,'smpd')

jun_edfd <- spatial_map(6,'June',2012,edfd_file,'eddi.D2')
jun_smpd <- spatial_map(6,'June',2012,smpd_file,'smpd')

jul_edfd <- spatial_map(7,'July',2012,edfd_file,'eddi.D2')
jul_smpd <- spatial_map(7,'July',2012,smpd_file,'smpd')

aug_edfd <- spatial_map(8,'August',2012,edfd_file,'eddi.D2')
aug_smpd <- spatial_map(8,'August',2012,smpd_file,'smpd')

# sep_edfd <- spatial_map(8,'September',2012,edfd_file,'eddi.D2')
# sep_smpd <- spatial_map(8,'September',2012,smpd_file,'smpd')


remove_title <- function(obj){
  
  out <- obj+
    ggtitle('')
  
  return(out)
}

#save legend
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#make the legend bigger
big_legend <-
aug_smpd+
guides(fill=guide_legend(reverse=T))+
theme(legend.key.size = unit(1.5, 'in'), #change legend key size
      legend.key.height = unit(4.2, 'in'), #change legend key height
      legend.key.width = unit(0.14, 'in'), #change legend key width
      legend.title = element_text(size=12), #change legend title font size
      legend.text = element_text(size=12,face='bold')) 

legend_ <- get_legend(big_legend)

plot1 <-
  ggarrange(apr_edfd,apr_smpd,
          may_edfd,may_smpd,
          jun_edfd,jun_smpd,
          jul_edfd,jul_smpd,
          aug_edfd,aug_smpd,
          nrow=5,ncol=2,legend = F,labels = c('A','F','B','G','C','H','D','I','E','J'))

#no title
# plot1 <-
#   ggarrange(remove_title(apr_edfd),remove_title(apr_smpd),
#             remove_title(may_edfd),remove_title(may_smpd),
#             remove_title(jun_edfd),remove_title(jun_smpd),
#             remove_title(jul_edfd),remove_title(jul_smpd),
#             remove_title(aug_edfd),remove_title(aug_smpd),
#             nrow=5,ncol=2,legend = F,labels = c('A','F','B','G','C','H','D','I','E','J'))


outplot_for_paper <-
  ggdraw()+
draw_plot(plot1, x=0.00,y=0,width=0.9,height=1)+
  draw_plot(legend_, x=0.905,y=0.011,width=0.1,height=0.987)

outplot_for_paper

#save plot
  ggsave(outplot_for_paper,filename = paste0(new_dir,'/Monthly_comparision_between_edfd_smpd_2012.tiff'),dpi = 300,
         height = 19.61, width = 11.3, units = 'in')
  

  # Spatial maps, case study 2017-COMPLETED -------------------------------------------
  
  
  edfd_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/edfd_month_average.csv')
  smpd_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/smpd_month_average.csv')
  
  #test
  # month=6
  # year=2012
  # file=edfd_file
  # month_name='June'
  # index_='eddi.D2'
  
  spatial_map <- function(month, month_name, year, file, index_){
    #find index of correct date
    index_num <- which(file$month==month & file$year==year)
    subset <- file[index_num,3:ncol(file)]
    
    #Location of csv files for eddi
    box_ <- paste0(lpath, index_) 
    setwd(box_)
    directory_create <- ('Spatial_distribution')
    dir.create(directory_create)
    
    #need 1 EDDI file to overlay new cluster information on top of
    outpath = paste(box_, directory_create, sep = '/')
    
    ascii_file <- ('19810807.asc')
    fh_ascii <- read.asciigrid(ascii_file)
    
    #Open mask
    mask_file <- mask_wrapper(index_ = index_)
    
    mask_fh <- read.csv(mask_file) %>% 
      dplyr::select(2)
    
    
    #Overlay cluster information on top of ascii file to not create a new map
    for (row in 1:ncol(subset)){
      if (mask_fh[row,1] == 0 || subset[row]==0){
        fh_ascii@data[row,1] = NA
      }else{
        fh_ascii@data[row,1] = subset[row]
      }
    }
    
    #Create a dataframe from a EDDI ascii file with clusters already added
    fh_raster <- ascii_wrapper(index_ = index_, fh_ascii = fh_ascii)
    
    colnames(fh_raster) = c("long", "lat", "Obs")
    
    #Map of the CONUS
    baseData <- map_data('state')

    #Set color pallette
    colfunc_map <- colorRampPalette(c('yellow','red'))
    #Set breaks
    breaks_spatial <<- seq(min(0), max(fh_raster$Obs, na.rm = T), by =0.25)
    
    
    
    # #Test for discrete color scale
    # spatial_discrete <- 
    #   
    #   ggplot()+
    #   geom_tile(data = fh_raster,aes(x=long, y=lat, fill = factor(Obs)))+
    #   theme(legend.title = element_blank(), legend.position = 'none')+
    #   geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    #   theme_bw()+
    #   coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    #   labs(title=title_name, x="Longitude", y="Latitude", fill = "Weeks in\n Flash Drought")+ 
    #   xlab(label = "Longitude")+
    #   theme(plot.title = element_text(hjust = 0.5, size=15)) +
    #   scale_fill_brewer(na.value='white', breaks=c(0,30,60,90,120,150,180,210,240,270,300),
    #                     palette = 'YlOrRd', )
    # 
    # 
    # 
    # 
    # +
    #   scale_fill_gradientn(colors = colfunc(100),
    #                        breaks = breaks_spatial,
    #                        labels= as.character(breaks_spatial), 
    #                        na.value = 'white')+
    #   theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))
    # 
    # 
    
    
    spatial_final <-
      ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=month_name, x="Longitude", y="Latitude", fill = "% of Time\n in Drought")+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=15)) +
      scale_fill_gradientn(colors = colfunc_map(40),
                           breaks = breaks_spatial,
                           labels= as.character(breaks_spatial), 
                           na.value = 'white')+
      theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))+
      scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
      scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
      theme(plot.title = element_text(size=16,face='bold'),
            axis.title.y =element_text(size=11))
      
    
    
    if (index_ == 'eddi.D2'){
      name_o <- 'EDFD'
    }else{
      name_o <- 'SMPD'
    }
    
    grob <- grobTree(textGrob(name_o, x=0.9,  y=0.05, hjust=0,
                              gp=gpar(col="black", fontsize=10)))
    
    outmap <- spatial_final+annotation_custom(grob)
    
    return(outmap)
  }
  
  # mar_edfd <- spatial_map(3,'March',2012,edfd_file,'eddi.D2')
  # mar_smpd <- spatial_map(3,'March',2012,smpd_file,'smpd')
  
  #Only do this analysis for reviewers comments
  # apr_edfd <- spatial_map(4,'April',2017,edfd_file,'eddi.D2')
  # apr_smpd <- spatial_map(4,'April',2017,smpd_file,'smpd')
  
  may_edfd <- spatial_map(5,'May',2017,edfd_file,'eddi.D2')
  may_smpd <- spatial_map(5,'May',2017,smpd_file,'smpd')
  
  jun_edfd <- spatial_map(6,'June',2017,edfd_file,'eddi.D2')
  jun_smpd <- spatial_map(6,'June',2017,smpd_file,'smpd')
  
  jul_edfd <- spatial_map(7,'July',2017,edfd_file,'eddi.D2')
  jul_smpd <- spatial_map(7,'July',2017,smpd_file,'smpd')
  
  aug_edfd <- spatial_map(8,'August',2017,edfd_file,'eddi.D2')
  aug_smpd <- spatial_map(8,'August',2017,smpd_file,'smpd')
  
  sep_edfd <- spatial_map(9,'September',2017,edfd_file,'eddi.D2')
  sep_smpd <- spatial_map(9,'September',2017,smpd_file,'smpd')
  
  
  remove_title <- function(obj){
    out <- obj+
      ggtitle('')
    
    return(out)
  }
  
  #save legend
  get_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  #make the legend bigger
  big_legend <-
    aug_smpd+
    guides(fill=guide_legend(reverse=T))+
    theme(legend.key.size = unit(1, 'in'), #change legend key size
          legend.key.height = unit(3.3, 'in'), #change legend key height
          legend.key.width = unit(0.14, 'in'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=12,face='bold')) 
  
  
  legend_ <- get_legend(big_legend)
  
  # plot1 <-
  #   ggarrange(may_edfd,may_smpd,
  #             jun_edfd,jun_smpd,
  #             jul_edfd,jul_smpd,
  #             aug_edfd,aug_smpd,
  #             nrow=4,ncol=2,legend = F,labels = c('A','E','B','F','C','G','D','H'))
  # 
  plot1 <-
    ggarrange(may_edfd,may_smpd,
              jun_edfd,jun_smpd,
              jul_edfd,jul_smpd,
              aug_edfd,aug_smpd,
              sep_edfd,sep_smpd,
              nrow=5,ncol=2,legend = F,labels = c('A','F','B','G','C','H','D','I','E'))
  #no title
  # plot1 <-
  #   ggarrange(remove_title(apr_edfd),remove_title(apr_smpd),
  #             remove_title(may_edfd),remove_title(may_smpd),
  #             remove_title(jun_edfd),remove_title(jun_smpd),
  #             remove_title(jul_edfd),remove_title(jul_smpd),
  #             remove_title(aug_edfd),remove_title(aug_smpd),
  #             nrow=5,ncol=2,legend = F,labels = c('A','F','B','G','C','H','D','I','E','J'))
  
  
  outplot_for_paper <-
    ggdraw()+
    draw_plot(plot1, x=0.00,y=0,width=0.9,height=1)+
    draw_plot(legend_, x=0.905,y=0.011,width=0.1,height=0.987)
  outplot_for_paper
  #save plot
  ggsave(outplot_for_paper,filename = paste0(new_dir,'/Monthly_comparision_between_edfd_smpd_2017.tiff'),dpi = 300,
         height = 15.4, width = 10.5, units = 'in')
  
  
  
# Number of Flash Drought Events lasting longer than 2 weeks- Spatial Distribution-COMPLETED -------------------------------
#Test for one week
# week_count <- 2
# column_name <- 'X2_week'
# index_ = 'smeddi.smpd'
#
week_count_function <- function(week_count, column_name, index_){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  #Read spatial file
  spatial2 <- read.csv(paste0(index_,'_Spatial_distribution_greater_than_n_weeks.csv'), header =T, row.names =1)
  columns <- colnames(spatial2)
  
  print(paste0('Working on ',column_name, ' for index ',index_,'.'))
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #Open mask
  mask_file <- mask_wrapper(index_)
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0) {
      fh_ascii@data[row,1] = NA
    }
    else{
      fh_ascii@data[row,1] = spatial2[row,column_name]
    }
  }
  
  #Map of the CONUS
  baseData <- map_data('state')
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_, fh_ascii)
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  title_name <- paste(week_count, ' weeks or greater', sep = '')
  
  #Keep track of all values
  if (str_detect(column_name, '2')){
    brek_week2 <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 3))
  }else if (str_detect(column_name, '3')){
    brek_week3 <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 3))
  }else if (str_detect(column_name, '4')){
    brek_week4 <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 3))
  }else if (str_detect(column_name, '5')){
    brek_week5 <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 3))
  }else if (str_detect(column_name, '6')){
    brek_week6 <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 2))
  }else if (str_detect(column_name, '7')){
    brek_week7 <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 1))
  }
  
  #Plot each graph, We only want to keep the limits from max of 2 week and min of 7 week
  if (str_detect(column_name, '2')){
    g2 <- ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_name, fill = "Number of Events")+ 
      xlab(label = "Longitude")+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5, size=15)) +
      scale_fill_gradientn(colors = colfunc_wyor(75),
                           breaks = seq(0,max(brek_week2),by =5),
                           labels= as.character(seq(0,max(brek_week2),by =5)), 
                           na.value = 'white',
                           limits = c(0,max(brek_week2))) +
      theme(legend.key.height = unit(2.0, 'cm'))
    
    return(g2)
    
  } else {
    g2 <- ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_name, fill = "Number of Events")+ 
      theme(plot.title = element_text(hjust = 0.5, size=15))+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank()) +
            scale_fill_gradientn(colors = colfunc_wyor(75),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(brek_week2))) +
      theme(legend.key.height = unit(0.8, 'cm'))
    return(g2)
  }
  
}

index_list <- c('eddi.D2','smpd')

# index_list <- c('smeddi.smpd')
for (index_ in index_list){
  out_path_weeks <- paste0(lpath, index_,'/Spatial_distribution')

   if (index_ == 'eddi.D2'){
    eddi.d2.s2<<- week_count_function(2,"X2_week",index_)
    eddi.d2.s3<<- week_count_function(3,"X3_week",index_)
    eddi.d2.s4<<- week_count_function(4,"X4_week",index_)
    eddi.d2.s5<<- week_count_function(5,"X5_week",index_)
    eddi.d2.s6<<- week_count_function(6,"X6_week",index_)
    eddi.d2.s7<<- week_count_function(7,"X7_week",index_)
    eddi.d2_final <<- ggarrange(eddi.d2.s2,eddi.d2.s3,
                               eddi.d2.s4,eddi.d2.s5,
                               eddi.d2.s6,eddi.d2.s7, ncol=3, nrow=2, common.legend = T, legend = c('right'))
    # #make a common title
    # annotate_figure(eddi.d2_final, top = text_grob(toupper(index_), face = 'bold', size = 24))
    # #save outputs
    # ggsave(filename = paste0(out_path_weeks,'/','event_count_distribution', '.tiff') , width = 10, height = 5,dpi = 300)
  } else if (index_ == 'smpd'){
    rzsm.s2<<- week_count_function(2,"X2_week",index_)
    rzsm.s3<<- week_count_function(3,"X3_week",index_)
    rzsm.s4<<- week_count_function(4,"X4_week",index_)
    rzsm.s5<<- week_count_function(5,"X5_week",index_)
    rzsm.s6<<- week_count_function(6,"X6_week",index_)
    rzsm.s7<<- week_count_function(7,"X7_week",index_)
    rzsm.final <<- ggarrange(rzsm.s2,rzsm.s3,rzsm.s4,
                            rzsm.s5,rzsm.s6,rzsm.s7,
                            ncol=3, nrow=2, common.legend = T, legend = c('right'))
    #make a common title
    annotate_figure(rzsm.final, top = text_grob(toupper(index_), face = 'bold', size = 24))
    #save outputs
    ggsave(filename = paste0(out_path_weeks,'/','event_count_distribution', '.tiff') , 
           width = 10, height = 5,dpi = 300)
  }   
}


#EDDI change scale of images
#put on same scale (max for eddi is 69 ... nice)
change_scale_eddi.d2_1 <- function(object, week_number){
  
  object_out <- object + labs(title = NULL) + 
    ylab(label=NULL)+ xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,69)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size=30, angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    theme(panel.border = element_blank())
  
  return(object_out)
}

change_scale_eddi.d2_2 <- function(object, week_number){
  
  object_out <- object + labs(title =NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,69)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    theme(panel.border = element_blank())
  
  return(object_out)
}

# change_scale_eddi.d2_1(eddi.d2.s2,2) +
#   theme(panel.border = element_blank())
# change_scale_eddi.d2_2(eddi.d2.s3,3)
# change_scale_eddi.d2_2(eddi.d2.s4,4)
# change_scale_eddi.d2_2(eddi.d2.s5,5)
# change_scale_eddi.d2_2(eddi.d2.s6,6)
# change_scale_eddi.d2_2(eddi.d2.s7,7)


#put on same scale (max for eddi is 69 ... nice)
change_scale_smpd_1.1 <- function(object){
  
  object_out <- object + labs(title=NULL) + ylab(label='SMPD')+ xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = seq(0,70,by=5),
                         labels= as.character(seq(0,70,by=5)),
                         na.value = 'white',
                         limits = c(0,69)) +
    theme(legend.key.height = unit(0.9, 'cm'),
          legend.key.width = unit(2.8, 'cm'),
          plot.title = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size=30, angle = 90),
          panel.border = element_blank(),
          legend.position="bottom",
          legend.spacing.x = unit(0, 'cm'),
          legend.margin=margin(),
          legend.title = element_text())+
          guides(fill = guide_legend(reverse = T,
                                     label.position = 'bottom',
                                     nrow=1,
                                     title = 'Number of \nEvents',title.vjust = 1.0))
  return(object_out)
}

#This will be only for a combined image for event count and event percentage
change_scale_final_image <- function(object){
  
  object_out <- object + labs(title=NULL) + ylab(label='SMPD')+ xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = seq(0,70,by=5),
                         labels= as.character(seq(0,70,by=5)),
                         na.value = 'white',
                         limits = c(0,69)) +
    theme(legend.key.height = unit(0.7, 'cm'),
          legend.key.width = unit(1.65, 'cm'),
          plot.title = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size=30, angle = 90),
          panel.border = element_blank(),
          legend.position="bottom",
          legend.spacing.x = unit(0, 'cm'),
          legend.margin=margin(),
          legend.title = element_text())+
    guides(fill = guide_legend(reverse = T,
                               label.position = 'bottom',
                               nrow=1,
                               title = 'Number of \nEvents',title.vjust = 1.0))
  return(object_out)
}
legend_final_image <- g_legend(change_scale_final_image(rzsm.s2))

count_legend <- g_legend(change_scale_smpd_1.1(rzsm.s2))

change_scale_smpd_1.2 <- function(object){
  
  object_out <- object + labs(title=NULL) + ylab(label=NULL)+ xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL,
                         na.value = 'white',
                         limits = c(0,69)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size=30, angle = 90))+
    theme(panel.border = element_blank())
  return(object_out)
}


change_scale_smpd_2 <- function(object){
  
  object_out <- object + labs(title=NULL) +
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,69)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_blank())+
    theme(panel.border = element_blank(),
          axis.ticks.y = element_blank())

  return(object_out)
}

#6 weeks code


combined_plots <- ggarrange(ncol=5,nrow =2,
          change_scale_eddi.d2_1(eddi.d2.s2,2),change_scale_eddi.d2_2(eddi.d2.s3,3),
          change_scale_eddi.d2_2(eddi.d2.s4,4),change_scale_eddi.d2_2(eddi.d2.s5,5),
          change_scale_eddi.d2_2(eddi.d2.s6,6),change_scale_smpd_1.2(rzsm.s2),
          change_scale_smpd_2(rzsm.s3),change_scale_smpd_2(rzsm.s4),change_scale_smpd_2(rzsm.s5),
          change_scale_smpd_2(rzsm.s6), font.label = list(size=20),vjust = 1.5)


event_week_count <- ggdraw() + 
  draw_plot(combined_plots, x=0.05,y=0.4,width=.95,height=0.5)+
  draw_plot(count_legend, x=.07,y=0.15,width=0.9,height=0.2)+
  draw_text("A", x = 0.02, y = 0.90,face='bold',size=30)+
  draw_text("EDFD", x = 0.040, y = 0.80,fontface='bold',size=18,angle=90)+
  draw_text("SMPD", x = 0.040, y = 0.55,fontface='bold',size=18,angle=90)+
  draw_text("N>2 weeks", x = 0.13, y = 0.909,face='bold',size=18)+
  draw_text("N>3 weeks", x = 0.32, y = 0.909,face='bold',size=18)+
  draw_text("N>4 weeks", x = 0.51, y = 0.909,face='bold',size=18)+
  draw_text("N>5 weeks", x = 0.70, y = 0.909,face='bold',size=18)+
  draw_text("N>6 weeks", x = 0.89, y = 0.909,face='bold',size=18)
  
event_week_count

ggsave(filename = paste0(new_dir,'/Event_count_by_WeekLength_edfd_smpd_no_winter.tiff'),dpi = 300,
       height = 7.61,width=16.7)

font_size_final = 14
event_week_count_final <- ggdraw() + 
  draw_plot(combined_plots, x=0.05,y=0.4,width=.95,height=0.5)+
  draw_plot(legend_final_image, x=.07,y=0.15,width=0.9,height=0.2)+
  draw_text("A", x = 0.02, y = 0.91,face='bold',size=23)+
  draw_text("EDFD", x = 0.055, y = 0.80,fontface='bold',size=font_size_final,angle=90)+
  draw_text("SMPD", x = 0.055, y = 0.55,fontface='bold',size=font_size_final,angle=90)+
  draw_text("N>2 weeks", x = 0.14, y = 0.909,face='bold',size=font_size_final)+
  draw_text("N>3 weeks", x = 0.33, y = 0.909,face='bold',size=font_size_final)+
  draw_text("N>4 weeks", x = 0.52, y = 0.909,face='bold',size=font_size_final)+
  draw_text("N>5 weeks", x = 0.71, y = 0.909,face='bold',size=font_size_final)+
  draw_text("N>6 weeks", x = 0.90, y = 0.909,face='bold',size=font_size_final)

# ggarrange(ncol = 1,nrow = 2,
#           combined_plots,count_legend,vjust = 4.0)
# 
# ggsave(filename = paste0(new_dir,'/Event_count_by_WeekLength_edfd_smpd.tiff'),dpi = 300,height = 6.63,width = 16.2) #Saving 10.3 x 4.56 in image
# 

# 
# #Old
# event_count_6_week_plot <- ggdraw()+
#   draw_plot(change_scale_eddi.d2_1(eddi.d2.s2,2), x = -.005, y = .4975, width = .2, height = .5 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s3,3), x = 0.2, y = .522, width = .2, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s4,4), x = 0.4, y = .522, width = .2, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s5,5), x = 0.60, y = .522, width = .2, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s6,6), x = 0.80, y = .522, width = .2, height = .45 ) +
#   draw_plot(change_scale_smpd_1.2(rzsm.s2), x = -.005, y = .022, width = .2, height = .43) +
#   draw_plot(change_scale_smpd_2(rzsm.s3), x = 0.2, y = .022, width = .2, height = .45 ) +
#   draw_plot(change_scale_smpd_2(rzsm.s4), x = 0.4, y = .022, width = .2, height = .45 ) +
#   draw_plot(change_scale_smpd_2(rzsm.s5), x = 0.6, y = .022, width = .2, height = .45 ) +
#   draw_plot(change_scale_smpd_2(rzsm.s6), x = 0.8, y = .022, width = .2, height = .45 ) +
#   draw_plot(count_legend,x=0,y=0)

#
# ggsave(filename = paste0(new_dir,'/Event_count_by_WeekLength_edfd_smpd.tiff'),dpi = 300) #Saving 10.3 x 4.56 in image

# ggsave(filename = paste0(new_dir,'/Event_count_by_WeekLength_edfd_smpd.tiff'),dpi = 300,plot = event_count_6_week_plot,
#        width =12.5, height=4.4,units='in')

#7 weeks code below
# ggdraw()+
#   draw_plot(change_scale_eddi.d2_1(eddi.d2.s2,2), x = -.005, y = .4975, width = .185, height = .5 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s3,3), x = 0.1727, y = .522, width = .17, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s4,4), x = 0.337, y = .522, width = .17, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s5,5), x = 0.5, y = .522, width = .17, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s6,6), x = 0.6635, y = .522, width = .17, height = .45 ) +
#   draw_plot(change_scale_eddi.d2_2(eddi.d2.s7,7), x = .826, y = .522, width = .17, height = .45 ) +
# draw_plot(change_scale_smpd_1.2(rzsm.s2), x = -.005, y = .022, width = .185, height = .45 ) +
# draw_plot(change_scale_smpd_2(rzsm.s3), x = 0.1727, y = .022, width = .17, height = .45 ) +
# draw_plot(change_scale_smpd_2(rzsm.s4), x = 0.337, y = .022, width = .17, height = .45 ) +
# draw_plot(change_scale_smpd_2(rzsm.s5), x = 0.5, y = .022, width = .17, height = .45 ) +
# draw_plot(change_scale_smpd_2(rzsm.s6), x = 0.6635, y = .022, width = .17, height = .45 ) +
# draw_plot(change_scale_smpd_2(rzsm.s7), x = .826, y = .022, width = .17, height = .45 ) +
#   draw_plot(count_legend,x=0,y=0)




# 
# ggarrange(change_scale_eddi.d2_1(eddi.d2.s2,2),
#           change_scale_eddi.d2_2(eddi.d2.s3,3),
#           change_scale_eddi.d2_2(eddi.d2.s4,4),
#           change_scale_eddi.d2_2(eddi.d2.s5,5),
#           change_scale_eddi.d2_2(eddi.d2.s6,6),
#           change_scale_eddi.d2_2(eddi.d2.s7,7),
#           change_scale_smpd_1.2(rzsm.s2),
#           change_scale_smpd_2(rzsm.s3),
#           change_scale_smpd_2(rzsm.s4),
#           change_scale_smpd_2(rzsm.s5),
#           change_scale_smpd_2(rzsm.s6),
#           change_scale_smpd_2(rzsm.s7),
#           common.legend = T,
#           legend = c('right'),
#           ncol=6,
#           nrow=2)

########## Seasonal Distribution Function for number of events that occurred in that season ##################### ---------------------------------------------------

#test column_name
column_name = 'MAM'
index_='eddi.D2'
#

list.files()

#Run function to find the seasonal distribution for each index
season_function <- function(column_name, index_){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  
  #Open mask
  mask_file <- mask_wrapper(index_)
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  #Read file for each index
  if (str_detect(index_,'smeddi')){
    file2 <- 'smeddi_seasonal_distribution.csv'
  }else if (index_== 'eddi.D2'){  
    file2 <- 'eddi.D2_seasonal_distribution.csv'
  }else{
    file2 <- paste0(index_,'_seasonal_distribution.csv')
  }
  
  fh <- read.csv(file2, row.names =1, header =T)
  fh_cols <- colnames(fh)
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0) {
      fh_ascii@data[row,1] = NA
    }else{
      fh_ascii@data[row,1] = fh[row,column_name]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_, fh_ascii)
  colnames(fh_raster) = c("long", "lat", "Obs")
  baseData <- map_data('state')
  
  #Find the maximum breaks
  if (column_name == 'MAM'){
    brek_MAM <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 8))
  }else if (column_name == 'JJA'){
    brek_JJA <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 8))
  }else if (column_name == 'SON'){
    brek_SON <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 8))
  }else if (column_name == 'DJF'){
    brek_DJF <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 8))
  }else if (column_name == 'MarNov'){
    brek_MARNOV <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 8))
  }
  
  if (column_name == 'MarNov'){
    title2 <- 'Mar-Nov'
  }else{
    title2 <- column_name
  }
  
  #Plot each graph
  colfunc <- colorRampPalette(c('yellow', 'red'))
  
  if (column_name == 'MarNov'){
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = "Total Weeks in\n Flash Drought")+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+ 
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(min(brek_MARNOV),max(brek_MARNOV))) +
      theme(legend.key.height = unit(1.2, 'cm'))
    
    return(g1) 
    
  }else if (column_name == 'JJA'){
    colfunc <- colorRampPalette(c('yellow', 'red'))
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = element_blank())+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = brek_JJA,
                           labels= as.character(brek_JJA), na.value = 'white') +
      theme(legend.key.height = unit(1.9, 'cm'))
    
    return(g1)
    
  } else {
    colfunc <- colorRampPalette(c('yellow', 'red'))
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = element_blank())+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(brek_JJA))) +
      theme(legend.key.height = unit(0.8, 'cm'))
    
    return(g1)
  }
}


##########Seasonal Distribution, when flash drought started by month - no month lag ###########

#test column_name
#jan, feb, march
column_name = 'JJA'
index_='eddi.D2'
#


#Run function to find the seasonal distribution for each index
season_month_function_no_lag <- function(index_, column_name){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  
  #Open mask
  mask_file <- mask_wrapper(index_)
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  mask_fh_t = t(mask_fh)
  
  #Read file for each index
  if (str_detect(index_,'smpd')){
    file2 <- 'smpd_FD_month_start_count.csv'
  }else if (index_== 'eddi.D2'){  
    file2 <- 'eddi.D2_FD_month_start_count.csv'
  }

  
  fh <- read.csv(file2, row.names =1, header =T)
  fh_cols <- colnames(fh)
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #overlay onto .asc file
  if (column_name == 'DJF'){
    range1=12
    range2=1
    range3=2
  } else if (column_name == 'MAM'){
    range_ = 3:5 
  } else if (column_name == 'JJA'){
    range_ = 6:8 
  } else if (column_name == 'SON'){
    range_ = 9:11 
  }

  #Can't write these 2 loops together, curly braces issues
  if (column_name !='DJF'){
  #overlay the sum of all months for each season
    for (grid in 1:ncol(mask_fh_t)){
      if (mask_fh_t[1,grid] == 0) {
        fh_ascii@data[grid,1] = NA
      } else {
      fh_ascii@data[grid,1] = sum(fh[range_,grid])
      }
    }
  }

  if (column_name == 'DJF'){
    for (grid in 1:ncol(mask_fh_t)){
      if (mask_fh_t[1,grid] == 0) {
        fh_ascii@data[grid,1] = NA
      }else{
        fh_ascii@data[grid,1] = fh[range1,grid] + fh[range2,grid] + fh[range3,grid]
      } 
    }
  }
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_, fh_ascii)
  colnames(fh_raster) = c("long", "lat", "Obs")
  baseData <- map_data('state')
  
  #Find the maximum breaks
  if (column_name == 'MAM'){
    brek_MAM <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 2))
  }else if (column_name == 'JJA'){
    brek_JJA <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 2))
  }else if (column_name == 'SON'){
    brek_SON <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 2))
  }else if (column_name == 'DJF'){
    brek_DJF <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 2))
  }
  #TITLE FOR PLOTS
  title2 <- column_name
  
  #Plot each graph
  colfunc <- colorRampPalette(c('yellow', 'red'))
  
  #highest plot will be summer season  (july, august, september)
  if (column_name == 'JJA'){
    colfunc <- colorRampPalette(c('yellow', 'red'))
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = element_blank())+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = brek_JJA,
                           labels= as.character(brek_JJA), na.value = 'white') +
      theme(legend.key.height = unit(1.9, 'cm'))
    
    return(g1)
    
  } else {
    colfunc <- colorRampPalette(c('yellow', 'red'))
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = element_blank())+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(brek_JJA))) +
      theme(legend.key.height = unit(0.8, 'cm'))
    
    return(g1)
  }
  }



# Call season month function - no month lag - COMPLETED ----------------------------------------------
index_list <- c('eddi.D2')
for (index_ in index_list){
  #create plots for each season
  jja2m <- season_month_function_no_lag(column_name = 'JJA',index_)
  mam2m <- season_month_function_no_lag(column_name = 'MAM',index_)
  son2m <- season_month_function_no_lag(column_name = 'SON',index_)
  # djf2m <- season_month_function_no_lag(column_name = 'DJF',index_)
}

# #Rounding was messing up the figure---OLD CODE WITH WINTER
# djf2m <- djf2m+
# scale_fill_gradientn(colors = colfunc_wyfr(100),
#                      breaks = NULL,
#                      labels= NULL, na.value = 'white',
#                      limits = c(0,max(max_value_season_eddi)+1))
#check the breaks to find scales, max is 80 for JJA
max_value_season_eddi <- max(brek_MAM, brek_JJA, brek_SON, na.rm=T)
min_value_season_eddi <- min(brek_MAM, brek_JJA, brek_SON, na.rm=T)

#jja2.eddi is the only legend to keep when graphing
jas_M2 <-  jja2m +
  scale_fill_gradientn(colors = colfunc_wyfr(100),
                       breaks = seq(0,max(max_value_season_eddi),by=2),
                       labels= as.character(seq(0,max(max_value_season_eddi),by=2)), na.value = 'white',
                       limits = c(0,max(max_value_season_eddi))) +
  theme(legend.key.height = unit(2.6, 'cm'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
        legend.text=element_text(size=10))+ 
  labs(fill='Number of Events \nby Season') 



#INDEX RZSM/SMPD Seasonal Distribution
index_list <- c('smpd')
for (index_ in index_list){
  #create plots for each season
  jja2Sm <- season_month_function_no_lag(column_name = 'JJA',index_)
  mam2Sm <- season_month_function_no_lag(column_name = 'MAM',index_)
  son2Sm <- season_month_function_no_lag(column_name = 'SON',index_)
  # djf2Sm <- season_month_function_no_lag(column_name = 'DJF',index_)
}

#change scale on outputs for everyting except JJA
#check the breaks to find scales, max is 80 for JJA
max_SM_season <- max(brek_MAM, brek_JJA, brek_SON, na.rm=T)
min_SM_season <- min(brek_MAM, brek_JJA, brek_SON, na.rm=T)


change_scale_smpd_MONTH <- function(object){
  object_out <- object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(min(min_SM_season),max(max_SM_season))) +
    theme(legend.key.height = unit(0.8, 'cm'),
          plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
          axis.ticks.x = element_blank(),
          panel.border = element_blank())
}

change_scale_smpd_jas <- function(object){
  object_out <- object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(min(min_SM_season),max(max_SM_season))) +
    theme(legend.key.height = unit(0.8, 'cm'),
          plot.title = element_text(hjust = 0.5, size=10, face = 'bold'),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

# djf_plot <- change_scale_smpd_MONTH(djf2Sm) + labs(title=NULL)
son_plot <- change_scale_smpd_MONTH(son2Sm)+ labs(title=NULL)
mam_plot <- change_scale_smpd_MONTH(mam2Sm)+ labs(title=NULL)
jja_plot <- change_scale_smpd_jas(jja2Sm)+ labs(title=NULL)


bold_text <- function(object){
  object_out <- 
    object +  theme(plot.title = element_text(face='bold'))
  
  return(object_out)
}

mam2m <-bold_text(mam2m)
jja2m <- bold_text(jas_M2)
son2m <- bold_text(son2m)
# djf2m <- bold_text(djf2m)

#Combine both EDDI and SMPD graphs
#Merge all plots together 
test <- ggarrange(mam2m,jja2m,son2m,
                  mam_plot,jja_plot, son_plot,
                  nrow = 2,  ncol=3,common.legend = T,  legend = c('right'))

ggdraw()+
  draw_plot(test,x=0.04,y=0,width = 0.95,height = 1)+
  draw_text("EDFD", x = 0.036, y = 0.78,fontface='bold',size=21,angle=90)+
  draw_text("SMPD", x = 0.036, y = 0.30,fontface='bold',size=21,angle=90)+
  draw_text("A", x = 0.02, y = 0.93,size=26)+
  draw_text("B", x = 0.02, y = 0.43,size=26)

ggsave(filename--  = paste0(new_dir,'/FD_by_month_edfd_smpd_no_lag.tiff'),
       dpi = 300,height=6.35,units='in',width=16.7)



# Find grid cell numbers --------------------------------------------------

#Location of csv files for eddi
box_ <- paste0(lpath, index_) 
setwd(box_)
#Open mask
if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
  mask_file <- ('NLDAS_mask_invert_No_puerto_rico_eddi.csv')
}else{
  mask_file <- ('Mask_final_puertorico_and_missing_values_rzsm.csv')
}

mask_fh <- read.csv(mask_file, row.names =1, header =T)
#need 1 EDDI file to overlay new cluster information on top of
fh_eddi <- read.asciigrid("19810807.asc")
head(fh_eddi)
a <- rev(fh_eddi$data)
a[59000]
fh_eddi$data[57101]
#Overlay cluster information on top of ascii file

file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/smpd/smpd_spatial_distribution.csv', header = T)
file_eddi <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/eddi.D2/eddi.D2_spatial_distribution.csv', header = T)

file[72000,2]
file_eddi[72000,2]
#EDDI.D2 index (find the grid pixels that I want)
for (row in 1:nrow(fh_eddi@data)){
    if(row > 23325 && row < 23369){
    fh_eddi@data[row,1] = file_eddi[row,2]
    } else {
    fh_eddi@data[row,1] = NA
  }
}

# #EDDI.D2 index (find the grid pixels that I want)
# for (row in 1:nrow(fh_eddi@data)){
#   if (row%%464 ==100 && row%%224==150){
#     fh_eddi@data[row,1] = row
#   } else {
#     fh_eddi@data[row,1] = NA
#   }
# }
# 
#   
#   
#   if (row > 23300 && row < 23350){
#     fh_eddi@data[row,1] = row
#   }  
#   if (row > 23300-464*2 && row < 23350-464*2){
#     fh_eddi@data[row,1] = row
#   }  
#   if (row > 23300-464*3 && row < 23350-464*3){
#     fh_eddi@data[row,1] = row
#   }  
#   if (row > 23300-464*4 && row < 23350-464*4){
#     fh_eddi@data[row,1] = row
#   }
#  
# 
# }


#Map of the CONUS
baseData <- map_data('state')

#Create a dataframe from a EDDI ascii file with clusters already added
if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
  fh_raster <- eddi_ras_wrap(fh_eddi)
}else{
  fh_raster <- else_ras_wrap(fh_eddi)
}

colnames(fh_raster) = c("long", "lat", "Obs")

colfunc2 <- colorRampPalette(c('white','yellow', 'red'))


#Plot each graph
 ggplot()+
  geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
  theme(legend.title = element_blank())+
  geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
  theme_bw()+
  coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
  labs(title=toupper(index_), x="Longitude", y="Latitude", fill = "% of Years with \n a Flash Drought")+ 
  theme(plot.title = element_text(hjust = 0.5, size=15),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank())+
   theme(legend.key.height = unit(1.5, 'cm'))






#potential Montana

##########  Percent of years with a flash drought -- COMPLETED ##########  -------------------------------------------------------
#function for graphs

percent_function <- function(index_){
  box_ <- paste0(lpath, index_) 
  setwd(box_)

  #First looking at all flash droughts
  if (index_ == 'eddi.D2'){
    csv_file <- paste0('percent_affected_eddi.csv')
  }else if (index_== 'smpd'){
    csv_file <- paste0('percent_affected_',index_,'.csv')
  }
  
  
  
  #Open mask
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    mask_file <- paste0(lpath,'NLDAS_mask_invert_No_puerto_rico_eddi.csv')
  }else{
    mask_file <- paste0(lpath,'Mask_final_puertorico_and_missing_values_rzsm.csv')
  }
  
  mask_fh <- read.csv(mask_file, row.names =1, header =T)
  fh <- vroom(csv_file) %>% 
    dplyr::select(., 2:ncol(.))
  
  #need 1 EDDI file to overlay new cluster information on top of
  fh_eddi <- read.asciigrid("19810807.asc")
  head(fh_eddi)
  
  #Overlay cluster information on top of ascii file
  for (row in 1:nrow(fh_eddi@data)){
    if (mask_fh[row,1] == 0) {
      fh_eddi@data[row,1] = NA
    }
    else{
      fh_eddi@data[row,1] = fh[row,1]
    }
  }
  
  #Map of the CONUS
  baseData <- map_data('state')
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    fh_raster <- eddi_ras_wrap(fh_eddi)
  }else{
    fh_raster <- else_ras_wrap(fh_eddi)
  }
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  colfunc2 <- colorRampPalette(c('white','yellow', 'red'))
  
  #Monitor the breaks for eddi defintions only
  
  #Keep track of all values
  if (index_ == 'eddi.D1'){
    brek_e1p <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 10))
  }else if (index_ == 'eddi.D2'){
    brek_e2p <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 10))
    max_edfd <<- max(fh_raster['Obs'], na.rm =T)
  }else if  (index_ == 'eddi.D3'){
    brek_e3p <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 10))
  }else if (index_ =='smpd'){
    max_smpd <<- max(fh_raster['Obs'], na.rm =T)
  }
  
  brek_percent <<- ceiling(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), by= 10))
  
  #Plot each graph
  percent_final <- ggplot()+
    geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
    theme(legend.title = element_blank())+
    geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    theme_bw()+
    coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    labs(title=toupper(index_), x="Longitude", y="Latitude", fill = "% of Years with \n a Flash Drought")+ 
    theme(plot.title = element_text(hjust = 0.5, size=15),
          axis.title = element_blank()) +
    scale_fill_gradientn(colors = colfunc(100),
                         breaks = brek_percent,
                         labels = as.character(brek_percent),
                         na.value = 'white',
                         limits = c(min(brek_percent),max(brek_percent)))+
    theme(legend.key.height = unit(1.5, 'cm'))+
    scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
    scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))
  
  percent_final
  
  return(percent_final)
}

ed2.per <- percent_function('eddi.D2')
rzsm.per <- percent_function('smpd')



#combine into 1 graph , each having its own scale
#max scale is smpd
max_smpd <- ceiling(max_smpd)

eddi_percent <- ed2.per+ scale_fill_gradientn(colors = colfunc_wyfr(300),
                                      breaks = NULL,
                                      labels = NULL,
                                      na.value = 'white',
                                      limits = c(min(0),max_smpd))+ labs(title='EDFD')+
  theme(legend.key.height = unit(1.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())




rzsm_percent <-
  rzsm.per+ scale_fill_gradientn(colors = colfunc_wyfr(100),
                               breaks = seq(0,max_smpd,by=10),
                               labels = as.character(seq(0,max_smpd,by=10)),
                               na.value = 'white',
                               limits = c(min(0),max_smpd))+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.text = element_text(size=10),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


#Combine
ggarrange(eddi_percent,rzsm_percent,common.legend = T,legend = c('right'),labels = c('A','B'),
          font.label = list(size=25,color='black'))


ggsave(filename = paste0(new_dir,'/Percent_of_years_edfd_smpd.tiff'),
       dpi = 300,height=5.54,width=16.7,units = 'in')


###########   Average fd length of occurence -- COMPLETED ########## -----------------------------------------------

avg_length_function <- function(index_){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  
  #overlay file
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #Open mask
  mask_file <- mask_wrapper(index_ = index_)
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    csv_file <- 'sum_of_occurences_eddi.csv'
  }else{
    csv_file <- paste0('sum_of_occurences_',index_,'.csv')
  }
  
  #open Sum_of_occurences file
  fh <- vroom(csv_file) %>% 
    dplyr::select(., 2:ncol(.))
  
  #Overlay information on top of ascii file
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0) {
      fh_ascii@data[row,1] = NA
    }
    else{
      fh_ascii@data[row,1] = fh[row,1]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    fh_raster <- eddi_ras_wrap(fh_ascii)
  }else{
    fh_raster <- else_ras_wrap(fh_ascii)
  }
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  #Keep track of all values of breaks (EDDI definitions)
  if (index_ == 'eddi.D1'){
    brek_e1sum <<- round(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), length.out= 8),digits = 2)
  }else if (index_ == 'eddi.D2'){
    brek_e2psum <<- round(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), length.out= 8),digits = 2)
  }else if  (index_ == 'eddi.D3'){
    brek_e3psum <<- round(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T), length.out= 8),digits = 2)
  } else if (index_ == 'smpd'){
    max_smpd <<-  max(fh_raster['Obs'], na.rm=T)
  }
  
  #Keep track of all values of breaks (all other definitions)
  brek_sum <- round(seq(min(fh_raster['Obs'], na.rm =T), max(fh_raster['Obs'], na.rm =T),length.out=8),digits = 2)
  
  #Map of the CONUS, base of ggplot
  baseData <- map_data('state')
  #set the color pallette for white =0, with also yellow and red
  colfunc2 <<- colorRampPalette(c('white','yellow', 'red'))
  #Plot each graph
  average_graph <- ggplot()+
    geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
    theme(legend.title = element_blank(), legend.position = 'none')+
    geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    theme_bw()+
    coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    labs(title=toupper(index_), x="Longitude", y="Latitude", fill = "Average Flash Drought \nWeeks per Year\n")+ 
    theme(plot.title = element_text(hjust = 0.5, size=15),
          axis.title = element_blank()) +
    scale_fill_gradientn(colors = colfunc2(60),
                         breaks = brek_sum,
                         labels = as.character(brek_sum), 
                         na.value = 'white') +
    theme(legend.key.height = unit(1.5, 'cm'))+
    scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
    scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))
  
  return(average_graph)
}

ed2.avg <- avg_length_function('eddi.D2') #maximum is 5.5
rzsm.avg <- avg_length_function('smpd') #maximum is 5.8


#combine into 1 graph (supplemental figure)
#max scale is 95 for smpd
eddi_average <- ed2.avg+ scale_fill_gradientn(colors = colfunc_wyfr(100),
                                              breaks = NULL,
                                              labels = NULL,
                                              na.value = 'white',
                                              limits = c(min(0),max_smpd))+ labs(title='EDFD')+
  theme(legend.key.height = unit(1.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())




rzsm_average <- rzsm.avg+ scale_fill_gradientn(colors = colfunc_wyfr(100),
                                               breaks = round(seq(0,max_smpd,length.out=8),digits = 2),
                                               labels = as.character(round(seq(0,max_smpd,length.out=8),digits = 2)),
                                               na.value = 'white',
                                               limits = c(min(0),max_smpd))+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.text = element_text(size=10),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


#Combine
ggarrange(eddi_average,rzsm_average,common.legend = T,legend = c('right'),labels = c('A','B'),
          font.label = list(size=25,color='black'))

ggsave(filename = paste0(new_dir,'/Average_weeks_by_year_edfd_smpd.tiff'),
       dpi = 300,height=5.54,width=16.7,units = 'in')
#Multi_magnitude
# Multi-magnitude occurence (Percent of Years with a flash drought) YR-----------------------------------------------
# #test
# category_number = 1
# index_ = 'rzsm'
# #
multi_graph_function <- function(category_number, index_, legend_spacing=9){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  
  #First looking at all flash droughts
  if (index_=='eddi.D2'){
    csv_file <- 'eddi.D2_Multi_magnitude_percent.csv'
  }else{
    csv_file <- paste0(index_,'_Multi_magnitude_percent.csv')
  }
  
  #index magnitude category file
  mult_file1 <- read.csv(csv_file,header=T, row.names=1)
  
  #Open mask
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    mask_file <- ('NLDAS_mask_invert_No_puerto_rico_eddi.csv')
  }else{
    mask_file <- ('Mask_final_puertorico_and_missing_values_rzsm.csv')
  }
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  fh_ascii <- read.asciigrid("19810807.asc")
  
  #Overlay cluster information on top of ascii file
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0) {
      fh_ascii@data[row,1] = NA
    }else{
      fh_ascii@data[row,1] = mult_file1[row,category_number]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    fh_raster <- eddi_ras_wrap(fh_ascii)
  }else{
    fh_raster <- else_ras_wrap(fh_ascii)
  }
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  title_plot <- paste0(toupper(index_),'_',as.character(category_number))
  
  #make a common legend for the final product with maximum a 
  if (category_number ==1 & index_=='eddi.D1'){
    ed1_mul_brk<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
    mul_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }else if (category_number ==1 & index_=='eddi.D2'){
    ed2_mul_brk<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
    mul_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }else if (category_number ==1 & index_=='eddi.D3'){
    ed3_mul_brk<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
    mul_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }else if (category_number ==1){
    mul_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }
  
  baseData <- map_data('state')
  colfunc2 <<- colorRampPalette(c('white','yellow', 'red'))
  
  #Plot each graph
  if (category_number == 1){
    multi_plot <- ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_plot, x="Longitude", y="Latitude", fill = "% of Years\nwith one F.D.")+ 
      
      theme(plot.title = element_text(hjust = 0.5, size=15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc2(100),
                           breaks = mul_brek,
                           labels= as.character(mul_brek),
                           na.value = 'white',limits = c(0,max(mul_brek))) +
      theme(legend.key.height = unit(1.0, 'in'),
            legend.key.width = unit(0.5, 'cm'))
    return(multi_plot)
    
  }else{
    multi_plot <- ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_plot, x="Longitude", y="Latitude", fill = "% of Years\nwith one F.D.")+ 
      theme(plot.title = element_text(hjust = 0.5, size=15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank())+
      scale_fill_gradientn(colors = colfunc2(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(mul_brek)))
    return(multi_plot)
  }
}

index_ ='eddi.D2'
#Call function Definition 2 EDDI
e2.ED1 <- multi_graph_function(category_number = 1, index_ = index_)
e2.ED2 <- multi_graph_function(category_number = 2, index_ = index_)
e2.ED3 <- multi_graph_function(category_number = 3, index_ = index_)
e2.ED4 <- multi_graph_function(category_number = 4, index_ = index_)

index_ ='smpd'
#Call function RZSM/SMPD
rzsm.ED1 <- multi_graph_function(category_number = 1, index_ = index_)
rzsm.ED2 <- multi_graph_function(category_number = 2, index_ = index_)
rzsm.ED3 <- multi_graph_function(category_number = 3, index_ = index_)
rzsm.ED4 <- multi_graph_function(category_number = 4, index_ = index_)

#max is 80 for eddi, but 90 for SM
eddi_d1 <- e2.ED1+ 
  scale_fill_gradientn(colors = colfunc_wyfr(100),
                       breaks = NULL,
                       labels = NULL,
                       na.value = 'white',
                       limits = c(min(0),max(90)))+ labs(title='ED1')+
  theme(legend.key.height = unit(1.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank())

#only for eddi category 2,3,4
eddi_mag_func <- function(object, category_number,index_){
 
  if (index_=='eddi'){
    cate <- 'ED'
  }else{
    cate <- 'SM'
  }
  
  object_out <-  object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels = NULL,
                         na.value = 'white',
                         limits = c(min(0),max(90)))+ 
   labs(title=paste0(cate, category_number))+
    theme(legend.key.height = unit(1.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
 return(object_out)
}

eddi_d2 <- eddi_mag_func(e2.ED2,2,'eddi')
eddi_d3 <- eddi_mag_func(e2.ED3,3,'eddi')
eddi_d4 <- eddi_mag_func(e2.ED4,4,'eddi')


smpd_d1 <- rzsm.ED1+ 
  scale_fill_gradientn(colors = colfunc_wyfr(100),
                       breaks = seq(0,90,by=5),
                       labels = as.character(seq(0,90,by=5)),
                       na.value = 'white',
                       limits = c(min(0),max(90)))+ labs(title='SM1')+
  theme(legend.key.height = unit(2.2, 'cm'),
        plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank())
  
          

smpd_d2 <- eddi_mag_func(rzsm.ED2,2,'smpd')
smpd_d3 <- eddi_mag_func(rzsm.ED3,3,'smpd')
smpd_d4 <- eddi_mag_func(rzsm.ED4,4,'smpd')


ggarrange(eddi_d1,eddi_d2,eddi_d3,eddi_d4,
          smpd_d1,smpd_d2,smpd_d3,smpd_d4,
          ncol=4,nrow=2,common.legend = T,
          legend = c('right'), labels = c('A','','','','B','','',''),
          font.label = list(size=20))

ggsave(filename = paste0(new_dir2,'/Percent_of_years_edfd_smpd.tiff'),
       dpi = 300,height=5.55,width=14.2,units='in')


# Multi-magnitude occurence (Percent of Weeks with a flash drought) WK-----------------------------------------------
# #test
# category_number = 1
# index_ = 'rzsm'
# #
multi_graph_function <- function(category_number, index_, legend_spacing=8){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  
  #First looking at all flash droughts
  if (index_=='eddi.D2'){
    csv_file <- 'eddi.D2_Multi_magnitude_weekly_percent.csv'
  }else{
    csv_file <- paste0(index_,'_Multi_magnitude_weekly_percent.csv')
  }
  
  #index magnitude category file
  mult_file1 <- read.csv(csv_file,header=T, row.names=1)
  
  #Open mask
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    mask_file <- ('NLDAS_mask_invert_No_puerto_rico_eddi.csv')
  }else{
    mask_file <- ('Mask_final_puertorico_and_missing_values_rzsm.csv')
  }
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  fh_ascii <- read.asciigrid("19810807.asc")
  
  #Overlay cluster information on top of ascii file
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0) {
      fh_ascii@data[row,1] = NA
    }else{
      fh_ascii@data[row,1] = mult_file1[row,category_number]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    fh_raster <- eddi_ras_wrap(fh_ascii)
  }else{
    fh_raster <- else_ras_wrap(fh_ascii)
  }
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  if (index_=='eddi.D2'){
    cat_name = 'ED'
  }else{
    cat_name='SM'
  }
  title_plot <- paste0(cat_name,as.character(category_number))
  
  #make a common legend for the final product with maximum a 
  if (category_number ==1 & index_=='eddi.D1'){
    ed1_mul_brk<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
    mul_brek<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
  }else if (category_number ==1 & index_=='eddi.D2'){
    ed2_mul_brk<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
    mul_brek<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
  }else if (category_number ==1 & index_=='eddi.D3'){
    ed3_mul_brk<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
    mul_brek<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
  }else if (category_number ==1){
    mul_brek<<- round(seq(0, max(fh_raster$Obs,na.rm=T), length.out=legend_spacing),digits=1)
  }
  
  baseData <- map_data('state')
  colfunc2 <<- colorRampPalette(c('white','yellow', 'red'))
  
  #Plot each graph
  if (category_number == 1){
    multi_plot <-
      ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_plot, x="Longitude", y="Latitude", fill = "% of Weeks\n in Flash Drought")+ 
      
      theme(plot.title = element_text(hjust = 0.5, size=15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc2(100),
                           breaks = mul_brek,
                           labels= as.character(mul_brek),
                           na.value = 'white',limits = c(0,max(mul_brek))) +
      theme(legend.key.height = unit(1.0, 'in'),
            legend.key.width = unit(0.5, 'cm'))
    return(multi_plot)
    
  }else{
    multi_plot <- ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_plot, x="Longitude", y="Latitude", fill = "% of Years\nwith one F.D.")+ 
      theme(plot.title = element_text(hjust = 0.5, size=15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank())+
      scale_fill_gradientn(colors = colfunc2(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(mul_brek)))
    return(multi_plot)
  }
}

index_ ='eddi.D2'
#Call function Definition 2 EDDI
e2.ED1 <-  multi_graph_function(category_number = 1, index_ = index_)
e2.ED2 <- multi_graph_function(category_number = 2, index_ = index_)
e2.ED3 <- multi_graph_function(category_number = 3, index_ = index_)
e2.ED4 <- multi_graph_function(category_number = 4, index_ = index_)

index_ ='smpd'
#Call function RZSM/SMPD
rzsm.ED1 <- multi_graph_function(category_number = 1, index_ = index_)
rzsm.ED2 <- multi_graph_function(category_number = 2, index_ = index_)
rzsm.ED3 <- multi_graph_function(category_number = 3, index_ = index_)
rzsm.ED4 <- multi_graph_function(category_number = 4, index_ = index_)


#SmPD is hightest, keep as the legend
smpd_d1 <-
  rzsm.ED1+
  scale_fill_gradientn(colors = colfunc_wyfr(100),
                       breaks = mul_brek,
                       labels= as.character(mul_brek),
                       na.value = 'white',limits = c(0,max(mul_brek)))+
  theme(legend.key.height = unit(2.2, 'cm'),
        plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size=10))


# 
# #max is 80 for eddi, but 90 for SM
# eddi_d1 <- e2.ED1+ 
#   scale_fill_gradientn(colors = colfunc_wyfr(100),
#                        breaks = NULL,
#                        labels = NULL,
#                        na.value = 'white',
#                        limits = c(min(0),max(90)))+ labs(title='ED1')+
#   theme(legend.key.height = unit(1.5, 'cm'),
#         plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.border = element_blank())

#only for eddi category 2,3,4
eddi_mag_func <- function(object, category_number,index_){
  
  if (index_=='eddi'){
    cate <- 'ED'
  }else{
    cate <- 'SM'
  }
  
  object_out <-  object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels = NULL,
                         na.value = 'white',
                         limits = c(min(0),7.4))+ 
    labs(title=paste0(cate, category_number))+
    theme(legend.key.height = unit(1.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
  return(object_out)
}

eddi_d1 <-eddi_mag_func(e2.ED1,1,'eddi')
eddi_d2 <- eddi_mag_func(e2.ED2,2,'eddi')
eddi_d3 <- eddi_mag_func(e2.ED3,3,'eddi')
eddi_d4 <- eddi_mag_func(e2.ED4,4,'eddi')


smpd_d2 <- eddi_mag_func(rzsm.ED2,2,'smpd')
smpd_d3 <- eddi_mag_func(rzsm.ED3,3,'smpd')
smpd_d4 <- eddi_mag_func(rzsm.ED4,4,'smpd')

# 
smpd
# 



ggarrange(eddi_d1,eddi_d2,eddi_d3,eddi_d4,
          smpd_d1,smpd_d2,smpd_d3,smpd_d4,
          ncol=4,nrow=2,common.legend = T,
          legend = c('right'), labels = c('A','','','','B','','',''),
          font.label = list(size=20))

ggsave(filename = paste0(new_dir2,'/Percent_of_weeks_multi_category_edfd_smpd.tiff'),
       dpi = 300,height=5.55,width=14.2,units='in')


# Sum of occurences multi-magnitude distribution graphs -------------------

###Function
# #TEST
# index_='eddi.D1'
# category_number=1
# #

sum_multi_func <- function(category_number, index_, mult_file1,legend_spacing=9){
  
  
  #Open mask
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    mask_file <- ('NLDAS_mask_invert_No_puerto_rico_eddi.csv')
  }else{
    mask_file <- ('Mask_final_puertorico_and_missing_values_rzsm.csv')
  }
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  fh_ascii <- read.asciigrid("19810807.asc")
  
  #Overlay cluster information on top of ascii file
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0) {
      fh_ascii@data[row,1] = NA
    }else{
      fh_ascii@data[row,1] = mult_file1[category_number,row]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  if ((str_detect(index_,'eddi')) & (index_ !='smeddi')){
    fh_raster <- eddi_ras_wrap(fh_ascii)
  }else{
    fh_raster <- else_ras_wrap(fh_ascii)
  }
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  baseData <- map_data('state')
  title_name <- paste0(toupper(index_),'_',as.character(category_number))
  
  colfunc2 <- colorRampPalette(c('white','yellow','red'))
  
  #make a common legend for the final product with maximum a 
  if (category_number ==1 & index_=='eddi.D1'){
    ed1_sum_brk<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
    sum_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }else if (category_number ==1 & index_=='eddi.D2'){
    ed2_sum_brk<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
    sum_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }else if (category_number ==1 & index_=='eddi.D3'){
    ed3_sum_brk<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
    sum_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }else if (category_number ==1){
    sum_brek<<- seq(0, max(fh_raster$Obs,na.rm=T), by=legend_spacing)
  }
  
  if (category_number == 1){
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_name, x="Longitude", y="Latitude", fill = "Total Weeks in\n Flash Drought")+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=15),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc2(100),
                           breaks = sum_brek,
                           labels= as.character(sum_brek), na.value = 'white',
                           limits = c(0,max(sum_brek))) +
      theme(legend.key.height = unit(1.0, 'in'), legend.key.width = unit(0.5, 'cm'))
    
    return(g1)
  }else{
    
    #Keep up with the breaks for each category of each index for scale
    if(category_number == 2){
      brek_2sum <<- ceiling(seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), by = 3))
    }else if(category_number ==3){
      brek_3sum <<- ceiling(seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), by = 3))
    }else if(category_number ==4){
      brek_4sum <<- ceiling(seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), by = 3))
    }
    
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_name, x="Longitude", y="Latitude", fill = "Total Weeks in\n Flash Drought")+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=15),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank()) +
      scale_fill_gradientn(colors = colfunc2(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(sum_brek))) +
      theme(legend.key.height = unit(0.2, 'in'), legend.key.width = unit(0.5, 'cm'))
    
    return(g1)
  }
}


# EDDI.D2 Multi-sum graphs ------------------------------------------------

index_ = 'eddi.D2'
#Location of csv files for eddi
box_ <- paste0(lpath, index_) 
setwd(box_)

#First looking at all flash droughts
if (index_=='eddi.D2'){
  csv_file <- 'eddi.D2_Multi_sum_of_occurences.csv'
}else{
  csv_file <- paste0(index_,'_Multi_sum_of_occurences.csv')
}

mult_file1 <- read.csv(csv_file,header=T, row.names=1)

sum.ED1.d2 <- sum_multi_func(category_number=1,index_ = index_, mult_file1=mult_file1) 
sum.ED2.d2 <- sum_multi_func(2,index_ = index_,mult_file1) 
sum.ED3.d2 <- sum_multi_func(3,index_ = index_,mult_file1) 
sum.ED4.d2 <- sum_multi_func(4,index_ = index_,mult_file1) 


#max is 100 for eddi, but 130 for SM
sum_eddi_d1 <- sum.ED1.d2+ 
  scale_fill_gradientn(colors = colfunc_wyfr(100),
                       breaks = NULL,
                       labels = NULL,
                       na.value = 'white',
                       limits = c(min(0),max(130)))+ labs(title='ED1')+
  theme(legend.key.height = unit(1.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank())

#only for eddi category 2,3,4
eddi_mag_func <- function(object, category_number,index_){
  
  if (index_=='eddi'){
    cate <- 'ED'
  }else{
    cate <- 'SM'
  }
  
  object_out <-  object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels = NULL,
                         na.value = 'white',
                         limits = c(min(0),max(90)))+ 
    labs(title=paste0(cate,category_number))+
    theme(legend.key.height = unit(1.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank())
  
  return(object_out)
}

sum_eddi_d2 <- eddi_mag_func(sum.ED2.d2,2,'eddi')
sum_eddi_d3 <- eddi_mag_func(sum.ED3.d2,3,'eddi')
sum_eddi_d4 <- eddi_mag_func(sum.ED4.d2,4,'eddi')



# RZSM/SMPD Multi-sum graphs ------------------------------------------------
index_ = 'smpd'
#Location of csv files for eddi
box_ <- paste0(lpath, index_) 
setwd(box_)
csv_file <- paste0(index_,'_Multi_sum_of_occurences.csv')

mult_file1 <- read.csv(csv_file,header=T, row.names=1)

sum.rzsm.1 <- sum_multi_func(1,index_ = index_,mult_file1) 
sum.rzsm.2 <- sum_multi_func(2,index_ = index_,mult_file1) 
sum.rzsm.3 <- sum_multi_func(3,index_ = index_,mult_file1) 
sum.rzsm.4 <- sum_multi_func(4,index_ = index_,mult_file1) 




#max is 100 for eddi, but 130 for SM
sum_SM1 <- sum.rzsm.1+ 
  scale_fill_gradientn(colors = colfunc_wyfr(130),
                       breaks = seq(0,130,by=10),
                       labels = as.character(seq(0,130,by=10)),
                       na.value = 'white',
                       limits = c(min(0),max(130)))+ labs(title='SM1')+
  theme(legend.key.height = unit(2.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size=26, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        legend.text=element_text(size=10))


sum_smpd_d2 <- eddi_mag_func(sum.rzsm.2,2,'smpd')
sum_smpd_d3 <- eddi_mag_func(sum.rzsm.3,3,'smpd')
sum_smpd_d4 <- eddi_mag_func(sum.rzsm.4,4,'smpd')


ggarrange(sum_eddi_d1, sum_eddi_d2, sum_eddi_d3, sum_eddi_d4,
          sum_SM1,sum_smpd_d2,sum_smpd_d3,sum_smpd_d4,
          common.legend = T, legend = 'right',
          ncol=4,nrow=2 ,labels = c('A','','','','B','','',''),
          font.label = list(size=20))

ggsave(filename = paste0(new_dir2,'/Sum_of_weeks_edfd_smpd.tiff'),
       dpi = 300,height=5.55,width=14.2,units='in')

########### Event Count Percentage Distribution -- COMPLETED ########### -------------------------------------

colfunc_wyfr <- colorRampPalette(c('white','yellow','firebrick1','darkred'))
#Read files
edfd <- read.csv(paste0(lpath,'eddi.D2/eddi.D2_Spatial_distribution_greater_than_n_weeks.csv'),
                 header=T,row.names = 1)

edfd <- edfd[,1:6]
smpd <- read.csv(paste0(lpath,'smpd/smpd_Spatial_distribution_greater_than_n_weeks.csv'),
                 header=T,row.names = 1)

#find percentage of flash droughts for each row
edfd_out <- edfd %>% 
  mutate(total = X2_week + X3_week+ X4_week + X5_week + X6_week + X7_week) %>% 
  mutate(wk2 = (X2_week/total)*100)%>% 
  mutate(wk3 = (X3_week/total)*100)%>% 
  mutate(wk4 = (X4_week/total)*100)%>% 
  mutate(wk5 = (X5_week/total)*100)%>% 
  mutate(wk6 = (X6_week/total)*100)%>% 
  mutate(wk7 = (X7_week/total)*100)

#find percentage of flash droughts for each row
smpd_out  <- smpd %>% 
  mutate(total = X2_week + X3_week+ X4_week + X5_week + X6_week + X7_week) %>% 
  mutate(wk2 = (X2_week/total)*100)%>% 
  mutate(wk3 = (X3_week/total)*100)%>% 
  mutate(wk4 = (X4_week/total)*100)%>% 
  mutate(wk5 = (X5_week/total)*100)%>% 
  mutate(wk6 = (X6_week/total)*100)%>% 
  mutate(wk7 = (X7_week/total)*100)

edfd_out <- edfd_out[,8:ncol(edfd_out)]
smpd_out <- smpd_out[,8:ncol(smpd_out)]

#Now plot
max(edfd_out, na.rm=T)
max(smpd_out, na.rm=T)

max_brek <<- 82.5

week_function <- function(index_, week_desired){
  
  if (index_ =='eddi'){
    fh = edfd_out
  }else{
    fh = smpd_out
  }
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #Open mask
  mask_file <- mask_wrapper(index_ = index_)
  
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  
  
  #Overlay cluster information on top of ascii file to not create a new map
  for (row in 1:nrow(fh_ascii@data)){
    if (mask_fh[row,1] == 0){
      fh_ascii@data[row,1] = NA
    }else{
      fh_ascii@data[row,1] = fh[row,week_desired-1]
    }
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_ = index_, fh_ascii = fh_ascii)
  
  colnames(fh_raster) = c("long", "lat", "percentage")
  
  #Map of the CONUS
  baseData <- map_data('state')
  #Set title of ggplot
  title_name <- paste('>=',week_desired,' weeks' , sep ='')
  
  #minimum for eddi.D1 is 10
  spatial_final <-  ggplot()+
    geom_tile(data = fh_raster,aes(x=long, y=lat, fill = percentage))+
    theme(legend.title = element_blank(), legend.position = 'none')+
    geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    theme_bw()+
    coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    labs(title=title_name, x="Longitude", y="Latitude", fill = "Percentage")+ 
    xlab(label = "Longitude")+
    theme(plot.title = element_text(hjust = 0.5, size=15),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank()) +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels= NULL,
                         limits = c(0,max_brek),
                         na.value = 'white')+
    theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))
  
  return(spatial_final)
  
}

#put on same scale, adjust for proper picture
change_scale_ed2<- function(object){
  
  object_out <- object + labs(title = NULL) + 
    ylab(label=NULL)+ xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,max_brek)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size=30, angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    theme(panel.border = element_blank())
  
  return(object_out)
}
change_scale_ed3thru6 <- function(object){
  
  object_out <- object + labs(title =NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,max_brek)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    theme(panel.border = element_blank())
  
  return(object_out)
}

#Call functions 
ed2 <- week_function(index_='eddi',week_desired = 2)
ed2_f <- change_scale_ed2(ed2)

ed3 <- week_function(index_='eddi',week_desired = 3)
ed3_f <- change_scale_ed3thru6(ed3)

ed4 <- week_function(index_='eddi',week_desired = 4)
ed4_f <- change_scale_ed3thru6(ed4)

ed5 <- week_function(index_='eddi',week_desired = 5)
ed5_f <- change_scale_ed3thru6(ed5)

ed6 <- week_function(index_='eddi',week_desired = 6)
ed6_f <- change_scale_ed3thru6(ed6)

#smpd editing
#put on same scale 

CHANGE_LEGEND_HEIGHT <- 0.7


change_scale_sm2 <- function(object){
  
  object_out <- object + labs(title=NULL) +xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = seq(0,max_brek,by=5),
                         labels= as.character(seq(0,max_brek,by=5)),
                         na.value = 'white',
                         limits = c(0,max_brek)) +
    theme(legend.key.height = unit(CHANGE_LEGEND_HEIGHT, 'cm'),
          legend.key.width = unit(1.30, 'cm'),
          plot.title = element_text(size = 12, face = "bold"),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          legend.position="bottom",
          legend.spacing.x = unit(0, 'cm'),
          legend.margin=margin(),
          legend.title = element_text())+
    guides(fill = guide_legend(reverse = T,
                               label.position = 'bottom',
                               nrow=1,
                               title = 'Percentage of \nEvents',title.vjust = 1.0))
  return(object_out)
}

sm2 <- week_function(index_='smpd',week_desired = 2)
sm2_l <- change_scale_sm2(sm2)
count_legend_percent <- g_legend(sm2_l)

change_scale_sm2.2 <- function(object){
  
  object_out <- object + labs(title=NULL) + ylab(label=NULL)+ xlab(label=NULL)+
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL,
                         na.value = 'white',
                         limits = c(0,max_brek)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size=30, angle = 90))+
    theme(panel.border = element_blank())
  return(object_out)
}

sm2_f <- change_scale_sm2.2(sm2_l)

change_scale_sm3thru6 <- function(object){
  
  object_out <- object + labs(title=NULL) +
    scale_fill_gradientn(colors = colfunc_wlsofd(70),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,max_brek)) +
    theme(legend.key.height = unit(0.8, 'cm'))+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_blank(),
          axis.title = element_blank())+
    theme(panel.border = element_blank(),
          axis.ticks.y = element_blank())
  
  return(object_out)
}

sm3 <- week_function(index_='smpd',week_desired = 3)
sm3_f <- change_scale_sm3thru6(sm3)

sm4 <- week_function(index_='smpd',week_desired = 4)
sm4_f <- change_scale_sm3thru6(sm4)

sm5 <- week_function(index_='smpd',week_desired = 5)
sm5_f <- change_scale_sm3thru6(sm5)

sm6 <- week_function(index_='smpd',week_desired = 6)
sm6_f <- change_scale_sm3thru6(sm6)


#6 weeks code
combined_plots_percentage <- ggarrange(ncol=5,nrow =2,
                           ed2_f,ed3_f,ed4_f,ed5_f,ed6_f,
                           sm2_f,sm3_f,sm4_f,sm5_f,sm6_f,
                           font.label = list(size=20),vjust = 1.5)


font_size=14


percentage_week_count <-  ggdraw() + 
  draw_plot(combined_plots_percentage, x=0.05,y=0.4,width=.97,height=0.47)+
  draw_plot(count_legend_percent, x=.07,y=0.15,width=0.9,height=0.2)+
  draw_text("B", x = 0.02, y = 0.91,face='bold',size=23)+
  draw_text("EDFD", x = 0.055, y = 0.775,fontface='bold',size=font_size,angle=90)+
  draw_text("SMPD", x = 0.055, y = 0.539,fontface='bold',size=font_size,angle=90)+
  draw_text("N=2 weeks", x = 0.140, y = 0.879,face='bold',size=font_size)+
  draw_text("N=3 weeks", x = 0.334, y = 0.879,face='bold',size=font_size)+
  draw_text("N=4 weeks", x = 0.525, y = 0.879,face='bold',size=font_size)+
  draw_text("N=5 weeks", x = 0.721, y = 0.879,face='bold',size=font_size)+
  draw_text("N=6 weeks", x = 0.915, y = 0.879,face='bold',size=font_size)
percentage_week_count

# ggsave(filename = paste0(new_dir,'/Event_percentage_edfd_smpd.tiff'),dpi = 300,height = 6.63,width = 16.2) #Saving 10.3 x 4.56 in image

###Combine 2 plots
test <- ggdraw() +
  draw_plot(event_week_count_final,x=0.0,y=0.5,height=0.5,width = 1)+
  draw_plot(percentage_week_count,x=0,y=0,height = 0.55,width=0.98)

test

ggsave(filename = paste0(new_dir,'/Combined_Event_&_percentage.tiff'),dpi = 300,height = 7.78,width = 10.3,plot = test) #Saving 10.3 x 4.56 in image

_# event_percentage_6_week_plot <- ggdraw()+
#   draw_plot(ed2_f, x = -.005, y = .4975, width = .2, height = .5 ) +
#   draw_plot(ed3_f, x = 0.2, y = .522, width = .2, height = .45 ) +
#   draw_plot(ed4_f, x = 0.4, y = .522, width = .2, height = .45 ) +
#   draw_plot(ed5_f, x = 0.60, y = .522, width = .2, height = .45 ) +
#   draw_plot(ed6_f, x = 0.80, y = .522, width = .2, height = .45 ) +
#   draw_plot(sm2_f, x = -.005, y = .022, width = .2, height = .43) +
#   draw_plot(sm3_f, x = 0.2, y = .022, width = .2, height = .45 ) +
#   draw_plot(sm4_f, x = 0.4, y = .022, width = .2, height = .45 ) +
#   draw_plot(sm5_f, x = 0.6, y = .022, width = .2, height = .45 ) +
#   draw_plot(sm6_f, x = 0.8, y = .022, width = .2, height = .45 ) +
#   draw_plot(count_legend_percent,x=0,y=0)





# EDFD value prior to SMPD drought----------------------------------------

outpath = '/home/kdl/Insync/OneDrive/eddi_and_sm/Final_graph_outputs/EDFD_lag_to_SMPD'


#lead are actually lags
eplus4_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/4_week_edfd_lead_from_smpd.csv')
eplus3_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/3_week_edfd_lead_from_smpd.csv')
eplus2_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_edfd_lead_from_smpd.csv')
eplus1_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_edfd_lead_from_smpd.csv')
eminus0_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/0_week_edfd_lag_from_smpd.csv')
eminus1_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_edfd_lag_from_smpd.csv')
eminus2_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_edfd_lag_from_smpd.csv')
eminus3_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/3_week_edfd_lag_from_smpd.csv')
eminus4_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/4_week_edfd_lag_from_smpd.csv')



#test
file=eplus1_80
lead=T
number=1
lag=F



eddi_comparison_map <- function(file,lead,lag,number){
  #find index of correct date
  
  data <- file$X0
  
  #Location of csv files for eddi

  
  #need 1 EDDI file to overlay new cluster information on top of
  outpath = ('/home/kdl/Insync/OneDrive/eddi_and_sm/Final_graph_outputs/EDFD_lag_to_SMPD')
  
  ascii_file <- ('/home/kdl/Insync/OneDrive/eddi_and_sm/19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #Open mask
  index_='smpd'
  mask_file <- mask_wrapper(index_ = index_)
  
  mask_fh <- read.csv(mask_file) %>%
    dplyr::select(2)
  
  
  #bin the files so that its a factor
  #Overlay cluster information on top of ascii file to not create a new map
  for (row in 1:length(data)){
    if (mask_fh[row,1] == 0 || data[row] < -100){
      fh_ascii@data[row,1] = data[row]
    }else{
        fh_ascii@data[row,1] = data[row]
    }
    }
  
  
  index_='smpd'
  # function(object){
  #
  #  test <-  raster(fh_ascii) %>%
  #    flip(., direction ='x') %>%
  #    as.data.frame(., xy=TRUE)
  # fh_raster <- test
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_ = index_, fh_ascii = fh_ascii)
  
  colnames(fh_raster) = c("long", "lat", "Obs")
  
  #Map of the CONUS
  baseData <- map_data('state')
  #Set color pallette
  colfunc <<- colorRampPalette(c('yellow', 'red'))
  #Set breaks
  breaks_spatial <- seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), length.out=6)
  # breaks_spatial <- seq(0.2, max(0.8, na.rm = T), length.out=4)
  
  # breaks_spatial <- seq(0.8,1.0, length.out)
  
  a <- fh_raster$Obs<0.2
  length(a[a==TRUE])
  
  if (lead==T){
    title_f = paste0(number,'-week EDFD lag')
  }else if (lag==T){
    title_f = paste0(number,'-week EDFD lead')
  }
  
  # spatial_final <-
    # ggplot()+
    # geom_tile(data = fh_raster,aes(x=long, y=lat, fill = (Obs)))+
    # scale_fill_gradient(low = "white", high = "green", na.value='white')+
    # theme(legend.title = element_blank(), legend.position = 'none')+
    # geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    # theme_bw()+
    # coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    # labs(title=title_f, x="Longitude", y="Latitude", fill = "%\nTrue Positive ")+
    # xlab(label = "Longitude")+
    # theme(plot.title = element_text(hjust = 0.5, size=15)) +
    # theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))+
    # scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
    # scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
    # theme(plot.title = element_text(size=16,face='bold'),
    #       axis.title.y =element_text(size=11))
    
    spatial_final <-
    ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = (Obs)))+
      scale_fill_gradient(low = "yellow", high = "green", na.value='white')+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title_f, x="Longitude", y="Latitude", fill = "True Positive (%) ")+
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=15)) +
      theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))+
      scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
      scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
      theme(plot.title = element_text(size=16,face='bold'),
            axis.title.y =element_text(size=11))+
      scale_fill_gradientn(colors = colfunc_wyfr(300),
                           breaks = breaks_spatial,
                           limits = c(0.1,.8),
                           labels= as.character(breaks_spatial),
                           
                           na.value = 'white')
  
    spatial_final
  return(spatial_final)
}

# p4_80 <-eddi_comparison_map(eplus4_80,lead=T,lag=F,number=4)
p3_80 <-eddi_comparison_map(eplus3_80,lead=T,lag=F,number=3)
p2_80 <-eddi_comparison_map(eplus2_80,lead=T,lag=F,number=2)
p1_80 <-eddi_comparison_map(eplus1_80,lead=T,lag=F,number=1)
m0_80 <-eddi_comparison_map(eminus0_80,lead=F,lag=T,number=0)
m1_80 <-eddi_comparison_map(eminus1_80,lead=F,lag=T,number=1)
m2_80 <-eddi_comparison_map(eminus2_80,lead=F,lag=T,number=2)
m3_80 <-eddi_comparison_map(eminus3_80,lead=F,lag=T,number=3)
# m4_80 <-eddi_comparison_map(eminus4_80,lead=F,lag=T,number=4)



#save legend
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


test_legend <- p1_80+
  guides(fill=guide_legend(title="True Positive\nRate (%)", reverse=T))+
  theme(legend.key.height = unit(1.9, 'in'), legend.key.width = unit(0.5, 'cm'),
        legend.text = element_text(size=12),
        legend.title=element_text(size=14))
test_legend
remove_legend <- function(obj){
  
  obj_out <- obj+
    theme(legend.position = 'none')
  return(obj_out)
}

legend_edfd <- get_legend(test_legend)

#remove legends
p4_80o <- remove_legend(p4_80)
p3_80o <-remove_legend(p3_80)
p2_80o <- remove_legend(p2_80)
p1_80o <-remove_legend(p1_80)
m0_80o <-remove_legend(m0_80)
m1_80o <-remove_legend(m1_80)
m2_80o <-remove_legend(m2_80)
m3_80o <-remove_legend(m3_80)
m4_80o <-remove_legend(m4_80)

# out_plot1 <- ggarrange(p4_80o,p3_80o,p2_80o,p1_80o,
#                            common.legend = F,
#                            ncol =4,
#                            labels = c('A','B','C','D'))

out_plot1 <- ggarrange(p3_80o,p2_80o,p1_80o,
                       common.legend = F,
                       ncol =3,
                       labels = c('A','B','C'))

out_plot2<- ggarrange(m1_80o,m2_80o,m3_80o,
                      common.legend = F,
                      ncol =3,
                      labels = c('E','F','G'))

out_plot3 <- ggarrange(m0_80o)



out_final <-
  ggdraw()+
  draw_plot(out_plot1, x = 0.0, y = .65, width = 0.9, height = 0.35)+
    draw_plot(out_plot3, x = 0.2, y = .35, width = 0.5, height = 0.3)+  
    draw_plot(out_plot2, x = 0.0, y = .0, width = 0.9, height = 0.35)+  
  draw_plot(legend_edfd, x = 0.95, y = 0, width = .02, height = 1)+
  draw_text('D', x=0.32,y=0.65,fontface='bold', size=15)
  

out_final

ggsave(out_final,filename = paste0(outpath,'/','EDFD_lead_lag_to_smpd.tiff'), 
       width = 16.5, height = 10.3,dpi = 300,units='in')




# EDDI value prior to SMPD drought ----------------------------------------

outpath = '/home/kdl/Insync/OneDrive/eddi_and_sm/Final_graph_outputs/EDDI_lag_to_SMPD'

eplus2_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_eddi_lead_from_smpd_top_80_percent.csv')
eplus1_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_eddi_lead_from_smpd_top_80_percent.csv')
eminus0_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/0_week_eddi_lag_from_smpd_top_80_percent.csv')
eminus1_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_eddi_lag_from_smpd_top_80_percent.csv')
eminus2_80 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_eddi_lag_from_smpd_top_80_percent.csv')

eplus2_50 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_eddi_lead_from_smpd_top_50_percent.csv')
eplus1_50 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_eddi_lead_from_smpd_top_50_percent.csv')
eminus0_50 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/0_week_eddi_lag_from_smpd_top_50_percent.csv')
eminus1_50 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_eddi_lag_from_smpd_top_50_percent.csv')
eminus2_50 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_eddi_lag_from_smpd_top_50_percent.csv')

eplus2_20 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_eddi_lead_from_smpd_top_20_percent.csv')
eplus1_20 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_eddi_lead_from_smpd_top_20_percent.csv')
eminus0_20 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/0_week_eddi_lag_from_smpd_top_20_percent.csv')
eminus1_20 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/1_week_eddi_lag_from_smpd_top_20_percent.csv')
eminus2_20 <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/2_week_eddi_lag_from_smpd_top_20_percent.csv')


eplus2_80$X0

#test
file=eplus2_80
percent=80
lead=T
number=2
lag=F

file=eplus2_80
eddi_comparison_map <- function(file,percent,lead,lag,number){
  #find index of correct date

  data <- file$X0

  #Location of csv files for eddi
  box_ <- paste0(lpath, index_)
  setwd(box_)
  directory_create <- ('Spatial_distribution')
  dir.create(directory_create)

  #need 1 EDDI file to overlay new cluster information on top of
  outpath = paste(box_, directory_create, sep = '/')

  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)

  #Open mask
  index_='smpd'
  mask_file <- mask_wrapper(index_ = index_)

  mask_fh <- read.csv(mask_file) %>%
    dplyr::select(2)


  #bin the files so that its a factor
  #Overlay cluster information on top of ascii file to not create a new map
  for (row in 1:length(data)){
    if (mask_fh[row,1] == 0 || data[row] < -100){
      fh_ascii@data[row,1] = NA
    }else{
      if (data[row] < 0.54){
      fh_ascii@data[row,1] = 0
      }else if (data[row] >=0.54 && (data[row] < 0.84)){
        fh_ascii@data[row,1] = 0.5 
      }else if (data[row] >=0.84 && (data[row] < 1.30)){
        fh_ascii@data[row,1] = 1.0 
      }else if (data[row] >=1.30 && (data[row] < 1.65)){
        fh_ascii@data[row,1] = 1.5 
      }else if (data[row] >=1.65 && data[row] < 2.06){
        fh_ascii@data[row,1] = 2.0 
      }else if (data[row] >=2.06){
        fh_ascii@data[row,1] = 2.5 
      }
    }
  }


  index_='smpd'
  # function(object){
  #
  #  test <-  raster(fh_ascii) %>%
  #    flip(., direction ='x') %>%
  #    as.data.frame(., xy=TRUE)
  # fh_raster <- test
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_ = index_, fh_ascii = fh_ascii)

  colnames(fh_raster) = c("long", "lat", "Obs")

  #Map of the CONUS
  baseData <- map_data('state')
  #Set color pallette
  colfunc <<- colorRampPalette(c('yellow', 'red'))
  #Set breaks
  breaks_spatial <- seq(min(fh_raster$Obs, na.rm = T), max(fh_raster$Obs, na.rm = T), length.out=6)


  if (lead==T){
    title_f = paste0(number,'-week EDDI lead')
  }else if (lag==T){
    title_f = paste0(number,'-week EDDI lag')
  }

  spatial_final <-
    ggplot()+
    geom_tile(data = fh_raster,aes(x=long, y=lat, fill = (Obs)))+
      scale_fill_gradient(low = "white", high = "darkred", na.value='white')+
    theme(legend.title = element_blank(), legend.position = 'none')+
    geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
    theme_bw()+
    coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
    labs(title=title_f, x="Longitude", y="Latitude", fill = "EDDI values")+
    xlab(label = "Longitude")+
    theme(plot.title = element_text(hjust = 0.5, size=15)) +
      theme(legend.key.height = unit(0.8, 'in'), legend.key.width = unit(0.5, 'cm'))+
      scale_y_continuous(breaks=c(30,35,40,45), labels=c('30N', '35N','40N','45N'))+
      scale_x_continuous(breaks=c(-120,-110,-100,-90,-80,-70), labels=c('120W','110W','100W','90W','80W','70W'))+
      theme(plot.title = element_text(size=16,face='bold'),
            axis.title.y =element_text(size=11))
    
  return(spatial_final)
}

p2_80 <-eddi_comparison_map(eplus2_80,80,lead=T,lag=F,number=2)
p1_80 <-eddi_comparison_map(eplus1_80,80,lead=T,lag=F,number=1)
m0_80 <-eddi_comparison_map(eminus0_80,80,lead=F,lag=T,number=0)
m1_80 <-eddi_comparison_map(eminus1_80,80,lead=F,lag=T,number=1)
m2_80 <-eddi_comparison_map(eminus2_80,80,lead=F,lag=T,number=2)

#save legend
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


test_legend <- p2_80+
theme(legend.key.height = unit(1.9, 'in'), legend.key.width = unit(0.5, 'cm'))

remove_legend <- function(obj){
  
  obj_out <- obj+
    theme(legend.position = 'none')
  return(obj_out)
}

legend_eddi <- get_legend(test_legend)

#remove legends
p2_80o <- remove_legend(p2_80)
p1_80o <-remove_legend(p1_80)
m0_80o <-remove_legend(m0_80)
m1_80o <-remove_legend(m1_80)
m2_80o <-remove_legend(m2_80)

out_file_plot <- ggarrange(p2_80o,p1_80o,m0_80o,m1_80o,m2_80o,
          common.legend = F,
          ncol =1,
          labels = c('A','B','C','D','E'))

out <- ggdraw()+
draw_plot(out_file_plot, x = 0.0, y = .00, width = 0.94, height = 1)+
  draw_plot(legend_eddi, x = 0.88, y = 0, width = .07, height = 1)

ggsave(out,filename = paste0(outpath,'/','20_percent_of_cases_eddi_to_smpd.tiff'), 
       width = 5.22, height = 12,dpi = 300,units='in')

#50%
p2_50 <-eddi_comparison_map(eplus2_50,50,lead=T,lag=F,number=2)
p1_50 <-eddi_comparison_map(eplus1_50,50,lead=T,lag=F,number=1)
m0_50 <-eddi_comparison_map(eminus0_50,50,lead=F,lag=T,number=0)
m1_50 <-eddi_comparison_map(eminus1_50,50,lead=F,lag=T,number=1)
m2_50 <-eddi_comparison_map(eminus2_50,50,lead=F,lag=T,number=2)

#save legend
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


test_legend <- p2_50+
  theme(legend.key.height = unit(1.9, 'in'), legend.key.width = unit(0.5, 'cm'))

remove_legend <- function(obj){
  
  obj_out <- obj+
    theme(legend.position = 'none')
  return(obj_out)
}

legend_eddi <- get_legend(test_legend)

#remove legends
p2_50o <- remove_legend(p2_50)
p1_50o <-remove_legend(p1_50)
m0_50o <-remove_legend(m0_50)
m1_50o <-remove_legend(m1_50)
m2_50o <-remove_legend(m2_50)

out_file_plot <- ggarrange(p2_50o,p1_50o,m0_50o,m1_50o,m2_50o,
                           common.legend = F,
                           ncol =1,
                           labels = c('A','B','C','D','E'))

out <- ggdraw()+
  draw_plot(out_file_plot, x = 0.0, y = .00, width = 0.94, height = 1)+
  draw_plot(legend_eddi, x = 0.88, y = 0, width = .07, height = 1)

ggsave(out,filename = paste0(outpath,'/','50_percent_of_cases_eddi_to_smpd.tiff'), 
       width = 5.22, height = 12,dpi = 300,units='in')


#20%
p2_20 <-eddi_comparison_map(eplus2_20,20,lead=T,lag=F,number=2)
p1_20 <-eddi_comparison_map(eplus1_20,20,lead=T,lag=F,number=1)
m0_20 <-eddi_comparison_map(eminus0_20,20,lead=F,lag=T,number=0)
m1_20 <-eddi_comparison_map(eminus1_20,20,lead=F,lag=T,number=1)
m2_20 <-eddi_comparison_map(eminus2_20,20,lead=F,lag=T,number=2)

#save legend
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


test_legend <- p2_20+
  theme(legend.key.height = unit(1.9, 'in'), legend.key.width = unit(0.5, 'cm'))

remove_legend <- function(obj){
  
  obj_out <- obj+
    theme(legend.position = 'none')
  return(obj_out)
}

legend_eddi <- get_legend(test_legend)

#remove legends
p2_20o <- remove_legend(p2_20)
p1_20o <-remove_legend(p1_20)
m0_20o <-remove_legend(m0_20)
m1_20o <-remove_legend(m1_20)
m2_20o <-remove_legend(m2_20)

out_file_plot <- ggarrange(p2_20o,p1_20o,m0_20o,m1_20o,m2_20o,
                           common.legend = F,
                           ncol =1,
                           labels = c('A','B','C','D','E'))

out <- ggdraw()+
  draw_plot(out_file_plot, x = 0.0, y = .00, width = 0.94, height = 1)+
  draw_plot(legend_eddi, x = 0.88, y = 0, width = .07, height = 1)

ggsave(out,filename = paste0(outpath,'/','80_percent_of_cases_eddi_to_smpd.tiff'), 
       width = 5.22, height = 12,dpi = 300,units='in')



#make a big plot
ggarrange(p2_80o,p1_80o,m0_80o,m1_80o,m2_80o,
p2_50o,p1_50o,m0_50o,m1_50o,m2_50o,
p2_20o,p1_20o,m0_20o,m1_20o,m2_20o,
ncol=5,nrow=5)

#make a big plot
big_plot <- ggarrange(m2_80o,m2_50o,m2_20o,
          m1_80o,m1_50o,m1_20o,
          m0_80o,m0_50o,m0_20o,
          p1_80o,p1_50o,p1_20o,
          p2_80o,p2_80o,p2_20o,
          ncol=3,nrow=5,
          labels = c('A (20%)','B (50%)','C (80%)','','','','','','','','','','','',''),
          font.label=list(size=12))



test_legend <-
  p2_20+
  theme(legend.key.height = unit(1.9, 'in'), legend.key.width = unit(0.5, 'cm'))+
    guides(fill=guide_legend(title="EDDI\nvalues", reverse=T))

legend_eddi <- get_legend(test_legend)


out <-
  ggdraw()+
  draw_plot(big_plot, x = 0.0, y = .00, width = 0.94, height = 1)+
  draw_plot(legend_eddi, x = 0.93, y = 0, width = .07, height = 1)

ggsave(out,filename = paste0(outpath,'/','80_50_20_percent_of_cases_eddi_to_smpd.tiff'), 
       width =10.9, height = 12,dpi = 300,units='in')



# New_stuff_not_done ------------------------------------------------------


regional_dir

mon_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/Regional_analysis_between_clusters/montana_time_series.csv')
cal_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/Regional_analysis_between_clusters/california_time_series.csv')
maine_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/Regional_analysis_between_clusters/maine_time_series.csv')
kan_file <- read.csv('/home/kdl/Insync/OneDrive/eddi_and_sm/Regional_analysis_between_clusters/kansas_time_series.csv')

par(mfrow=c(2,1))
plot(mon_file)
lines(mon_file)
title('Montana - EDDI')
mtext(paste('Pearson R',cor(mon_file$X0, cal_file$X0, method=c('pearson')),'MK=',MannKendall(mon_file$X0)$tau[1]), sep=' ')
plot(cal_file)
lines(cal_file)
title('California - EDDI')
mtext(paste('MK=',MannKendall(cal_file$X0)$tau[1]),sep=' ')

savePlot()

cor(mon_file$X0, cal_file$X0, method=c('pearson'))

plot(maine_file)

MannKendall(mon_file$X0)$tau[1]

MannKendall(as.numeric(test[i,1:nrow(mon_file)]))


# # Seasonal Distribution --- OLD DON'T RUN  -------------------------------------------------
# index_list <- c('eddi.D2')
# for (index_ in index_list){
#   #create plots for each season
#   jja2 <- season_function(column_name = 'JJA',index_)
#   mam2 <- season_function(column_name = 'MAM',index_)
#   son2 <- season_function(column_name = 'SON',index_)
#   # djf2 <- season_function(column_name = 'DJF',index_)
#   # marnov2 <- season_function('MarNov',index_)
# }
# 
# #check the breaks to find scales from season_function
# max_value_season_eddi <- max( brek_MAM, brek_JJA, brek_SON, na.rm=T)
# min_value_season_eddi <- min( brek_MAM, brek_JJA, brek_SON, na.rm=T)
# 
# #highest value of 104 in JJA, rounding in season function messed it up
# max_value_season_eddi <- 104
# 
# 
# 
# #change scale on outputs
# change_scale_eddi.d2 <- function(object){
#   object_out <- object +
#     scale_fill_gradientn(colors = colfunc_wyfr(100),
#                          breaks = NULL,
#                          labels= NULL, na.value = 'white',
#                          limits = c(0,max(max_value_season_eddi))) +
#     theme(legend.key.height = unit(0.8, 'cm'),
#           axis.ticks.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           panel.border = element_blank(),
#           plot.title = element_text(hjust = 0.5, size=20, face = 'bold'))
#   
# }
# 
# #change scale on outputs
# change_scale_eddi.d2_MAM_only <- function(object){
#   object_out <- object +
#     scale_fill_gradientn(colors = colfunc_wyfr(100),
#                          breaks = NULL,
#                          labels= NULL, na.value = 'white',
#                          limits = c(0,max(max_value_season_eddi))) +
#     theme(legend.key.height = unit(0.8, 'cm'),
#           plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
#           axis.ticks.x = element_blank(),
#           axis.text.x = element_blank(),
#           panel.border = element_blank())
# }
# 
# #jja2.eddi is the only legend
# jja2.eddi <- jja2 +
#   scale_fill_gradientn(colors = colfunc_wyfr(100),
#                        breaks = seq(min_value_season_eddi,max(max_value_season_eddi),by=8),
#                        labels= as.character(seq(min_value_season_eddi,max(max_value_season_eddi),by=8)), na.value = 'white',
#                        limits = c(min_value_season_eddi,max(max_value_season_eddi))) +
#   theme(legend.key.height = unit(2.6, 'cm'),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         panel.border = element_blank(),
#         plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
#         legend.text=element_text(size=10))+ 
#     labs(fill='Number of Weeks\nin Flash Drought') 
# 
# # djf2.eddi <- change_scale_eddi.d2(djf2)
# mam2.eddi <- change_scale_eddi.d2_MAM_only(mam2)
# son2.eddi <- change_scale_eddi.d2(son2)
# 
# ggarrange(mam2.eddi,jja2.eddi,son2---OLD DONT RUN.eddi)
# 
# #INDEX RZSM/SMPD Seasonal Distribution
# index_list <- c('smpd')
# for (index_ in index_list){
#   #create plots for each season
#   jja.rzsm <- season_function('JJA',index_)
#   mam.rzsm <- season_function('MAM',index_)
#   son.rzsm <- season_function('SON',index_)
#   djf.rzsm <- season_function('DJF',index_)
#   marnov.rzsm <- season_function('MarNov',index_)
# }
# #check the breaks to find scales, max is 80 for JJA
# max_value_season <- max(brek_DJF, brek_MAM, brek_JJA, brek_SON, na.rm=T)
# min_value_season <- min(brek_DJF, brek_MAM, brek_JJA, brek_SON, na.rm=T)
# 
# #change scale on outputs for everyting except JJA
# 
# change_scale_smpd_mam_only <- function(object){
#   object_out <- object +
#     scale_fill_gradientn(colors = colfunc_wyfr(100),
#                          breaks = NULL,
#                          labels= NULL, na.value = 'white',
#                          limits = c(min(min_value_season),max(max_value_season))) +
#     theme(legend.key.height = unit(0.8, 'cm'),
#           plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
#           axis.ticks.x = element_blank(),
#           panel.border = element_blank())
# }
# change_scale_smpd_no_jja <- function(object){
#   object_out <- object +
#     scale_fill_gradientn(colors = colfunc_wyfr(100),
#                          breaks = NULL,
#                          labels= NULL, na.value = 'white',
#                          limits = c(min(min_value_season),max(max_value_season))) +
#     theme(legend.key.height = unit(0.8, 'cm'),
#           plot.title = element_text(hjust = 0.5, size=10, face = 'bold'),
#           panel.border = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_blank())
# }
# 
# jja.rzsm.a <- change_scale_smpd_no_jja(jja.rzsm) + labs(title=NULL)
# djf.rzsm.a <- change_scale_smpd_no_jja(djf.rzsm)+ labs(title=NULL)
# son.rzsm.a <- change_scale_smpd_no_jja(son.rzsm)+ labs(title=NULL)
# mam.rzsm.a <- change_scale_smpd_mam_only(mam.rzsm)+ labs(title=NULL)
# 
# 
# #Combine both EDDI and SMPD graphs
# #Merge all plots together 
# test <- ggarrange(mam2.eddi,jja2.eddi,son2.eddi,djf2.eddi,
#           mam.rzsm.a,jja.rzsm.a, son.rzsm.a,djf.rzsm.a,
#   nrow = 2,  ncol=4,common.legend = T,  legend = c('right'))
# 
# ggdraw()+
#   draw_plot(test,x=0.04,y=0,width = 0.95,height = 1)+
#   draw_text("EDFD", x = 0.036, y = 0.78,fontface='bold',size=21,angle=90)+
#   draw_text("SMPD", x = 0.036, y = 0.30,fontface='bold',size=21,angle=90)+
#   draw_text("A", x = 0.02, y = 0.93,size=26)+
#   draw_text("B", x = 0.02, y = 0.43,size=26)
#   
# ggsave(filename = paste0(new_dir,'/Seasonal_dist_edfd_smpd.tiff'),
#        dpi = 300,height=6.35,units='in',width=16.7)
# 
# # ggsave(filename = paste0(new_dir,'/Seasonal_dist_edfd_smpd.tiff'),
# #        dpi = 300)
# # 

##########Seasonal Distribution, when flash drought started by month - 1month lag ###########

#test column_name
#jan, feb, march
column_name = 'JAS'
index_='smpd'
#


#Run function to find the seasonal distribution for each index
season_month_function <- function(index_, column_name){
  
  #Location of csv files for eddi
  box_ <- paste0(lpath, index_) 
  setwd(box_)
  
  #Open mask
  mask_file <- mask_wrapper(index_)
  mask_fh <- read.csv(mask_file) %>% 
    dplyr::select(2)
  mask_fh_t = t(mask_fh)
  
  #Read file for each index
  if (str_detect(index_,'smpd')){
    file2 <- 'smpd_FD_month_start_count.csv'
  }else if (index_== 'eddi.D2'){  
    file2 <- 'eddi.D2_FD_month_start_count.csv'
  }
  
  fh <- read.csv(file2, row.names =1, header =T)
  fh_cols <- colnames(fh)
  
  ascii_file <- ('19810807.asc')
  fh_ascii <- read.asciigrid(ascii_file)
  
  #overlay onto .asc file
  if (column_name == 'JFM'){
    range_ = 1:3 
  } else if (column_name == 'AMJ'){
    range_ = 4:6 
  } else if (column_name == 'JAS'){
    range_ = 7:9 
  } else if (column_name == 'OND'){
    range_ = 10:12 
  }
  
  #overlay the sum of all months for each season
  for (grid in 1:ncol(mask_fh_t)){
    if (mask_fh_t[1,grid] == 0) {
      fh_ascii@data[grid,1] = NA
    }else{
      fh_ascii@data[grid,1] = sum(fh[range_,grid])
    } 
  }
  
  #Create a dataframe from a EDDI ascii file with clusters already added
  fh_raster <- ascii_wrapper(index_, fh_ascii)
  colnames(fh_raster) = c("long", "lat", "Obs")
  baseData <- map_data('state')
  
  #Find the maximum breaks
  if (column_name == 'JFM'){
    brek_JFM <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 5))
  }else if (column_name == 'AMJ'){
    brek_AMJ <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 5))
  }else if (column_name == 'JAS'){
    brek_JAS <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 5))
  }else if (column_name == 'OND'){
    brek_OND <<- ceiling(seq(min(fh_raster[,'Obs'],na.rm=T),max(fh_raster[,'Obs'], na.rm = T), by = 5))
  }
  
  #TITLE FOR PLOTS
  title2 <- column_name
  
  #Plot each graph
  colfunc <- colorRampPalette(c('yellow', 'red'))
  
  #highest plot will be summer season  (july, august, september)
  if (column_name == 'JAS'){
    colfunc <- colorRampPalette(c('yellow', 'red'))
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = element_blank())+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = brek_JAS,
                           labels= as.character(brek_JAS), na.value = 'white') +
      theme(legend.key.height = unit(1.9, 'cm'))
    
    return(g1)
    
  } else {
    colfunc <- colorRampPalette(c('yellow', 'red'))
    g1 <-  ggplot()+
      geom_tile(data = fh_raster,aes(x=long, y=lat, fill = Obs))+
      theme(legend.title = element_blank(), legend.position = 'none')+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)+
      theme_bw()+
      coord_fixed(ratio = 1.3,xlim = range(baseData$long), range(ylim = baseData$lat), expand = F, clip = 'on')+
      labs(title=title2, x="Longitude", y="Latitude", fill = element_blank())+ 
      xlab(label = "Longitude")+
      theme(plot.title = element_text(hjust = 0.5, size=18)) +
      scale_fill_discrete(na.value = 'white')+
      theme(axis.title.x = element_blank(), 
            axis.title.y  = element_blank(),
            axis.text=element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())+
      scale_fill_gradientn(colors = colfunc(100),
                           breaks = NULL,
                           labels= NULL, na.value = 'white',
                           limits = c(0,max(brek_JAS))) +
      theme(legend.key.height = unit(0.8, 'cm'))
    
    return(g1)
  }
}


# Call season month function ----------------------------------------------
index_list <- c('eddi.D2')
for (index_ in index_list){
  #create plots for each season
  jas_M <- season_month_function(column_name = 'JAS',index_)
  ond_M <- season_month_function(column_name = 'OND',index_)
  jfm_M <- season_month_function(column_name = 'JFM',index_)
  amj_M <- season_month_function(column_name = 'AMJ',index_)
}

#check the breaks to find scales from season_function
max_value_season_eddi <- max(brek_OND, brek_AMJ, brek_JAS, brek_JFM, na.rm=T)
min_value_season_eddi <- min(brek_OND, brek_AMJ, brek_JAS, brek_JFM, na.rm=T)


#jja2.eddi is the only legend
jas_M2 <- jas_M +
  scale_fill_gradientn(colors = colfunc_wyfr(100),
                       breaks = seq(0,max(max_value_season_eddi),by=2),
                       labels= as.character(seq(0,max(max_value_season_eddi),by=2)), na.value = 'white',
                       limits = c(0,max(max_value_season_eddi))) +
  theme(legend.key.height = unit(2.6, 'cm'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
        legend.text=element_text(size=10))+ 
  labs(fill='Number of Events \nby Season') +
  ggtitle('Jul-Sept')


#change scale on outputs
change_scale_month <- function(object){
  
  if (object$labels$title == 'OND'){
    out_title <- 'Oct-Dec'
  } else if (object$labels$title == 'JFM'){
    out_title <- 'Jan-Mar'
  } else if (object$labels$title == 'AMJ'){
    out_title <- 'Apr-Jun'
  }
  
  object_out <- object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(0,max(max_value_season_eddi))) +
    theme(legend.key.height = unit(0.8, 'cm'),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size=20, face = 'bold'))+
    ggtitle(out_title)
  
}

ond_M2 <- change_scale_month(ond_M)
jfm_M2 <- change_scale_month(jfm_M)
amj_M2 <- change_scale_month(amj_M)


ggarrange(amj_M2,jas_M2,ond_M2,jfm_M2, common.legend = T, legend = 'right')

#INDEX RZSM/SMPD Seasonal Distribution
index_list <- c('smpd')
for (index_ in index_list){
  #create plots for each season
  jas_SM <- season_month_function(column_name = 'JAS',index_)
  ond_SM <- season_month_function(column_name = 'OND',index_)
  jfm_SM <- season_month_function(column_name = 'JFM',index_)
  amj_SM <- season_month_function(column_name = 'AMJ',index_)
}

#check the breaks to find scales, max is 80 for JJA
max_SM_season <- max(brek_OND, brek_AMJ, brek_JAS, brek_JFM, na.rm=T)
min_SM_season <- min(brek_OND, brek_AMJ, brek_JAS, brek_JFM, na.rm=T)

#change scale on outputs for everyting except JJA

change_scale_smpd_MONTH <- function(object){
  object_out <- object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(min(min_SM_season),max(max_SM_season))) +
    theme(legend.key.height = unit(0.8, 'cm'),
          plot.title = element_text(hjust = 0.5, size=20, face = 'bold'),
          axis.ticks.x = element_blank(),
          panel.border = element_blank())
}

change_scale_smpd_jas <- function(object){
  object_out <- object +
    scale_fill_gradientn(colors = colfunc_wyfr(100),
                         breaks = NULL,
                         labels= NULL, na.value = 'white',
                         limits = c(min(min_SM_season),max(max_SM_season))) +
    theme(legend.key.height = unit(0.8, 'cm'),
          plot.title = element_text(hjust = 0.5, size=10, face = 'bold'),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

amj_plot <- change_scale_smpd_MONTH(amj_SM) + labs(title=NULL)
jfm_plot <- change_scale_smpd_MONTH(jfm_SM)+ labs(title=NULL)
ond_plot <- change_scale_smpd_MONTH(ond_SM)+ labs(title=NULL)
jas_plot <- change_scale_smpd_jas(jas_SM)+ labs(title=NULL)



ggarrange(amj_M2,jas_M2,ond_M2,jfm_M2, common.legend = T, legend = 'right')


#Combine both EDDI and SMPD graphs
#Merge all plots together 
test <- ggarrange(amj_M2,jas_M2,ond_M2,jfm_M2,
                  amj_plot,jas_plot, ond_plot,jfm_plot,
                  nrow = 2,  ncol=4,common.legend = T,  legend = c('right'))

ggdraw()+
  draw_plot(test,x=0.04,y=0,width = 0.95,height = 1)+
  draw_text("EDFD", x = 0.036, y = 0.78,fontface='bold',size=21,angle=90)+
  draw_text("SMPD", x = 0.036, y = 0.30,fontface='bold',size=21,angle=90)+
  draw_text("A", x = 0.02, y = 0.93,size=26)+
  draw_text("B", x = 0.02, y = 0.43,size=26)

ggsave(filename = paste0(new_dir,'/FD_by_month_edfd_smpd.tiff'),
       dpi = 300,height=6.35,units='in',width=16.7)

