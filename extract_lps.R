# code written by Paul Tanger, Brook Moyers or John Lovell. Please cite appropriately.

extract.lps<-function(model.in, stat, phe){
  toext<-model.in[[phe]]
  all.lps<-attr(toext, "lodprofile")
  lps.df<-data.frame()
  for(j in 1:length(all.lps)){
    qtl.id.in<-stat$qtl.id[stat$phenotype==phe][j]
    lp.df<-as.data.frame(all.lps[j])
    colnames(lp.df)<-c("chr", "pos", "lod")
    lp.df$phenotype<-phe
    lp.df$qtl.id<-qtl.id.in
    lp.df$st.lod<-lp.df$lod/max(lp.df$lod)
    lp.df$lod[lp.df$lod<0]<-NA
    lp.df$st.lod[lp.df$st.lod<0]<-NA
    lps.df<-rbind(lps.df,lp.df)
  }
  return(lps.df)
}