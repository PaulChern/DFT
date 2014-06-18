# density<-read.table( file = '/home/aju/Development/DFT/test/density_matrix.data')
# density<-as.matrix( density ) 
# 
# An<-read.table( file = '/home/aju/Development/DFT/test/enhancement.data')
# An<-apply( An, c(1,2), FUN = as.numeric )
# An[is.na(An)]<-0
# 
# 
# T<-read.table( file = '/home/aju/Development/DFT/test/kinetic.data' )
# T<-as.matrix( T )
# 
# K<-read.table( file = '/home/aju/Development/DFT/test/distance.data' )
# K<-as.matrix( K )
# 
# U<-read.table( file = '/home/aju/Development/DFT/test/deschol.data' )
# U<-as.matrix( U )

r<-read.table( file = '/home/aju/Development/DFT/test/radial_grid.data' )
p<-read.table( file = '/home/aju/Development/DFT/test/radial_density.data' )
plot( r[,1], p[,1], cex = 0.5, col = 'purple4' )

# image( density, col = rainbow(100) )
# image( An, col = rainbow(100) )
# image( rho, col = heat.colors(100) )
# 
# max(An)
# min(An)
# 
# max(K)
# min(K)
# 
# max(U)
# min(U)

#___________________________________________________________________________________________________
# library(rgl)
# rho<-read.table( file = '/home/aju/Development/DFT/test/rho.data' )
# rho<-as.matrix( rho )
# 
# r<-read.table( file = '/home/aju/Development/DFT/test/radial_grid.data' )
# r<-as.matrix( r )
# 
# th<-read.table( file = '/home/aju/Development/DFT/test/angular_grid.data' )
# th<-as.matrix( th )
# 
# d<-dim(rho)
# R<-matrix( r[,1], d[1], d[2], byrow = TRUE )
# TH<-matrix( th[,1], d[1], d[2], byrow = TRUE )
# PH<-matrix( seq( 0, 2*pi, length.out = d[1] ), d[1], d[2] )
# COL<-matrix( c( rep('gold', 40 * d[1] ), 
#                 rep('dodgerblue3', 10 * d[1] ),
#                 rep('black', 10 * d[1] ) ), d[1], d[2] )
# 
# 
# x<-rho * cos(TH) * cos(PH)
# y<-rho * cos(TH) * sin(PH)
# z<-rho * sin(TH)
# 
# 
# open3d()
# persp3d( r[,1], th[,1], rho, color = 'gold', alpha = 0.5, theta = 60 )
# 
# open3d()
# persp3d( x, y, z, color = COL, alpha = 0.5, add = TRUE )
# 
# 
# df<-data.frame(x=runif(1000,0,1),
#                y=runif(1000,0,1),
#                z=runif(1000,0,1),
#                color=round(runif(1000,1,10)))
# plot3d(df$x, df$y, df$z, col=df$color, size=2, type='s')
