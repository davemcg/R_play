ggplot(data=gt,aes(x=Pos,y=CCGO59heterozygosity)) + facet_wrap(~chrom,ncol=2) + 
  geom_segment(data=CCGO59,aes(x=start,xend=end,y="CCGO59",yend="CCGO59"),size=2,colour='red') +
  geom_segment(data=CCGO60,aes(x=start,xend=end,y="CCGO60",yend="CCGO60"),size=2,colour='orange') +
  geom_segment(data=CCGO61,aes(x=start,xend=end,y="CCGO61",yend="CCGO61"),size=2,colour='green') +
  geom_segment(data=CCGO62,aes(x=start,xend=end,y="CCGO62",yend="CCGO62"),size=2,colour='blue') +
  geom_segment(data=CCGO159,aes(x=start,xend=end,y="CCGO159",yend="CCGO159"),size=2,colour='purple') +
  geom_segment(data=CCGO160,aes(x=start,xend=end,y="CCGO160",yend="CCGO160"),size=2,colour='gray') +
  geom_segment(data=CCGO161,aes(x=start,xend=end,y="CCGO161",yend="CCGO161"),size=2,colour='turquoise')

