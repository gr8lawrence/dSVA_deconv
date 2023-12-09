require(MuSiC)
my_music2_prop_t_statistics = function(bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, samples, select.ct, expr_low=20, 
                                    prop_r=0.1, eps_c=0.05, eps_r=0.01, n_resample=20, sample_prop=0.5,cutoff_expr=0.05, 
                                    cutoff_fc=2, cutoff_c=0.05, cutoff_r=0.01, maxiter = 200, markers = NULL, 
                                    cell_size = NULL, ct.cov = FALSE, centered = FALSE, normalize = FALSE){
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if(length(gene.bulk) < 0.1*min(nrow(bulk.control.mtx), nrow(bulk.case.mtx))){
    stop('Not enough genes for bulk data! Please check gene annotations.')
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, ])
  
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if(length(gene_all) < 0.2*min(length(gene.bulk), nrow(sc.sce))){
    stop('Not enough genes between bulk and single-cell data! Please check gene annotations.')
  }
  
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  
  # remove lowly expressed genes from DE analysis: i.e., gene with average expression < expr_low
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr>=expr_low])
  
  # Analyse separately based on their disease status
  
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)];
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)];
  
  # Step 1: cell type deconvolution, set initial value
  # estimate cell type proportion for controls using music.
  
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce,
                            clusters = clusters, samples = samples, select.ct = select.ct,
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                            nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx =bulk.case, sc.sce = sc.sce,
                             clusters = clusters, samples = samples, select.ct = select.ct,
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                             nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control, prop_CASE)
  
  # start iteration
  iter=1
  ncell = length(select.ct)
  id_conv = NULL
  while(iter <= maxiter){
    # print(iter)
    # step 2: identify cell-type-specific DE genes
    # calculate mean/ sd of log fold change as an indicator of differential expression
    LOGFC = NULL
    MOD0=MOD1=matrix(0L, nrow = ncell, ncol = length(exp_genel))
    
    for(i in 1:n_resample){
      id_h = sample(colnames(bulk.control),round(ncol(bulk.control)*sample_prop))
      control_s = bulk.control[exp_genel, colnames(bulk.control) %in% id_h]
      prop_h = prop_control[colnames(control_s),]
      
      mod0 = apply(control_s, 1, function(x){
        mod = nnls(prop_h,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD0=MOD0+mod0
      
      id_d = sample(colnames(bulk.case),round(ncol(bulk.case)*sample_prop))
      case_s = bulk.case[exp_genel, colnames(bulk.case) %in% id_d]
      prop_d = prop_CASE[colnames(case_s),]
      
      mod1 = apply(case_s, 1, function(x){
        mod = nnls(prop_d,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD1=MOD1+mod1
      LOGFC = rbind(LOGFC,log1p(mod1)-log1p(mod0))
    }
    
    rcv=NULL
    mfc=NULL
    for(i in 1:ncell){
      s = LOGFC[seq(from=i,to=nrow(LOGFC), by=ncell),]
      rcv = rbind(rcv,apply(s,2,function(x){ifelse(mean(x)==0, 0, mean(x)/sd(x))}))
      mfc = rbind(mfc,apply(s,2,function(x){mean(x)}))
    }
    abs_rcv_logfc = abs(rcv)
    MOD0 = MOD0/n_resample
    MOD1 = MOD1/n_resample
    rownames(MOD0)=rownames(MOD1)=rownames(abs_rcv_logfc)=rownames(mfc)=select.ct
    
    # average cell type proportion
    mex = apply(prop_all,2,mean)
    lr = NULL
    for(celltype in select.ct){
      m = mex[celltype]
      rh = MOD0[celltype,]
      rd = MOD1[celltype,]
      fc = mfc[celltype,]
      # for genes with average expression within lower cutoff_expr for both conditions are removed from cell-type-specific DE genes detection
      llr = unique(intersect(names(rd[rd <= quantile(rd,prob=cutoff_expr)]),
                             names(rh[rh <= quantile(rh,prob=cutoff_expr)])))
      x = abs_rcv_logfc[celltype,]
      x = x[!names(x) %in% llr]
      # select genes with large mean/SD of log fc and with at least moderate log fc as DE
      if(m >= prop_r){
        lr = c(lr, intersect(names(x[x >= quantile(x,prob=1-cutoff_c)]), names(fc[abs(fc)>=log1p(cutoff_fc)])))
      }else{
        lr = c(lr, intersect(names(x[x >= quantile(x,prob=1-cutoff_r)]), names(fc[abs(fc)>=log1p(cutoff_fc)])))
      }
    }
    lr = unique(lr)
    
    # step 3: update sc gene list
    # remove identified DE genes from sc rna-seq data
    l = setdiff(gene_all,lr)
    sc.iter.sce = sc.sce[l,]
    
    # step 1: update cell type proportion based on new gene list
    if(length(id_conv)>0){
      case_sample = bulk.case[ , !colnames(bulk.case) %in% id_conv]
    }else{
      case_sample = bulk.case
    }
    
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce,
                           clusters=clusters, samples=samples, select.ct=select.ct,
                           markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                           nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize,verbose = F)$Est.prop.weighted
    
    prop_CASE = rbind(prop_case, prop_case_fix)
    
    if(length(id_conv)==1){
      rownames(prop_CASE) = c(rownames(prop_case), id_conv)
    }
    
    prop_all = rbind(prop_control,prop_CASE)
    
    # check convergence, by cell type
    prop_case=prop_case[rownames(prop_case_ini),]
    pc = abs(prop_case-prop_case_ini)
    conv = pc
    conv[,] = 1
    # use difference if rare cell type
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    # use percent change if common cell type
    pc[prop_case_ini > prop_r] = pc[prop_case_ini>prop_r]/prop_case_ini[prop_case_ini>prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > prop_r] < eps_c,0,1)
    convf = apply(conv,1,function(x){all(x==0)})
    
    # if an id converged, not updating anymore
    all_converge=FALSE
    id_conv = c(id_conv,names(convf[convf==TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% id_conv,]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv,]
    
    # if all converged or if only one subjects not converging--> music2 converged
    if(is.vector(prop_case_ini)){
      all_converge=TRUE
      break
    }else if(nrow(prop_case_ini)==0){
      all_converge=TRUE
      break
    }
    iter=iter+1
  }
  # return estimated proportion
  if(all_converge){
    return(list('Est.prop' = prop_all,'convergence'=TRUE,'n.iter'=iter,'DE.genes'=lr))}
  else{
    return(list('Est.prop' = prop_all,'convergence'=FALSE,'id.not.converge'=rownames(prop_case_ini)))}
}
