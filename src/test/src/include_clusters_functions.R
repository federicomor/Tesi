cat(crayon::cyan("\nGenerated these functions:\n"))

library(grid)
library(gridBase)
library(gridExtra)

################################################################################
#### mode correction
mode_correct_clusters = function(cl_old, cl_cur,verbose=0,very_verbose=0){ # and returns cl_new the updated cl_cur
	if(is.null(cl_old)){
		warning("The cl_old vector is NULL (probably it's just iteration 1 initialization).\nReturning original clustering.\n")
		return(cl_cur) # nothing to correct
	}
	if(very_verbose==1){
		verbose=1
	}
	
	cl_new = cl_cur
	new_labels = c()
	tied_labels = c() # solve them later
	already_taken_labels = c() # solve them later
	
	################################
	### Main loop ##################
	################################
	# for (i in 1:max(cl_cur)){
	for (i in unique(cl_cur)){
		indici_label_i = which(cl_cur == i)
		freq = table(cl_old[indici_label_i])
		if(very_verbose==1){
			cat("Cur label i =",i,"has mode table (^=index, _=occurences)")
			print(freq)
		}
		better_labels = as.numeric(names(freq[freq == max(freq)]))
		num_mode = length(better_labels)
		
		if(num_mode>=2){ # per risolvere i pareggi
			# warning("Tie case happened.\n")
			tied_labels = c(tied_labels,i)
			if(very_verbose==1){
				cat(crayon::magenta("Tied label case. Soving it later.\n"))
			}
			
		} else { # una sola moda, ma c'è da vedere se non è già presa
			if( !(better_labels %in% new_labels) ){
				better_label = better_labels
				new_labels = c(new_labels,better_label)
				cl_new[indici_label_i] = better_label
				if(very_verbose==1){
					cat("Assigning cur label",i,"to new label",better_label,"\n")
				}
			} else {
				already_taken_labels = c(already_taken_labels,i)
				if(very_verbose==1){
					cat(crayon::magenta("Already taken label case. Soving it later.\n"))
				}
			}
		}
		# cat(new_labels,"\n")
	}
	if(verbose==1){
		cat(crayon::green("##############\n"))
		cat("tied_labels =",tied_labels,"\n")
		cat("already_taken_labels =",already_taken_labels,"\n")
		cat("Available labels =",setdiff(unique(cl_cur),new_labels),"\n")
		cat(crayon::green("##############\n"))
	}
	
	
	################################
	### Tied labels ################
	################################
	if(very_verbose==1){
		cat(crayon::cyan("Solving tied labels\n"))
	}
	something_changed = 1
	while(length(tied_labels)!=0 && something_changed==1){
		something_changed = 0
		for(i in tied_labels){
			indici_label_i = which(cl_cur == i)
			freq = table(cl_old[indici_label_i])
			better_labels = as.numeric(names(freq[freq == max(freq)]))
			
			better_labels = setdiff(better_labels,new_labels)
			if(length(better_labels)>=1){
				# risolvi il pareggio ed eilimina label i da tied_labels
				better_label = better_labels[1]
				new_labels = c(new_labels,better_label)
				cl_new[indici_label_i] = better_label
				tied_labels = setdiff(tied_labels,i)
				if(very_verbose==1){
					cat("Assigning tied-label",i,"to new label",better_label,"\n")
					something_changed=1
				}
			}
			if(length(better_labels)==0){
				if(very_verbose==1){
					cat("Already taken label case. Soving it later.\n")
				}
				already_taken_labels = c(already_taken_labels,i)
			}
		}
	}
	
	
	################################
	### Already taken labels #######
	################################
	if(very_verbose==1){
		cat(crayon::cyan("Solving already taken labels\n"))
	}
	for(i in already_taken_labels){
		indici_label_i = which(cl_cur == i)
		# cat(indici_label_i,"\n")
		# for(k in 1:length(unique(cl_cur))){
		for(k in 1:max(cl_cur)){
			# cat(k,"\n")
			if( !(k %in% new_labels) ){
				new_labels = c(new_labels,k)
				cl_new[indici_label_i] = k
				if(very_verbose==1){
					cat("Assigning already-taken-label",i,"to new label",k,"\n")
				}
				break
			}
		}
	}
	if(verbose==1){
		cat("Cur clusters values = ",unique(cl_cur),"\n")
		cat("New clusters values = ",unique(cl_new),"\n")
		cat("cl_old =",cl_old,"\n")
		cat("cl_cur =",cl_cur,"\n")
		cat("cl_new =",cl_new,"\n")
	}
	if (length(unique(cl_new))!=length(unique(cl_new)) ||
		length(levels(factor(cl_new)))!=length(levels(factor(cl_cur))) ){
		cat(crayon::red("Something went wrong. Returning original clustering.\n"))
		return(cl_cur)
	}

	if(!all(cl_new==cl_cur)){
		cat(crayon::italic("Some change was made!\n"))
	}
	return(cl_new)
	
}
cat(crayon::red("- mode_correct_clusters(cl_old, cl_cur)\n"))

# cat(crayon::blue("Test sets\n"))
# ### test sets
# 
# cat(crayon::blue("- Ettore test\n"))
# cl_old = c(1,2,1,1,2,2,3,3)
# cl_cur = c(1,1,2,2,1,1,3,3)
# check_ = c(2,2,1,1,2,2,3,3)
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# # mode_correct_clusters(cl_old,cl_cur) == check_
# 
# cat(crayon::blue("- No change case\n"))
# cl_old = c(1,2,1,1,2,2,3,3)
# cl_cur = c(1,2,1,1,2,2,3,3)
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# 
# 
# cat(crayon::blue("- Tie case\n"))
# # tie case
# #   units  1 2 3 4 5 6 7 8 9 10
# cl_old = c(1,2,1,1,2,1,3,3)
# cl_cur = c(1,1,2,2,1,1,3,3)
# # ora (pareggio) vedere sia il cluster 1 è entrato in 2 ma anche viceversa
# check1 = c(2,2,1,1,2,2,3,3)
# check2 = c(1,1,2,2,1,1,3,3)
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# # mode_correct_clusters(cl_old,cl_cur) == check1
# # mode_correct_clusters(cl_old,cl_cur) == check2
# 
# cat(crayon::blue("- Double tie case\n"))
# # tie case double
# cl_old = c(1,2,1,1,2,1,3,3,4,4,5,5)
# cl_cur = c(1,1,2,2,1,1,3,3,4,5,4,5)
# # ora (pareggio) vedere sia il cluster 1 è entrato in 2 ma anche viceversa
# # e idem per cluster 4 e 5. Vediamo cosa ritorna
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# 
# cat(crayon::blue("- disappeared clusters test\n"))
# # case of collapsed clusters
# cl_old = c(1,1,1,2,2,2,1,1,3,4,4,4,4)
# cl_cur = c(1,1,1,1,1,1,1,1,1,2,2,2,2) 
# mode_correct_clusters(cl_old,cl_cur,verbose=1) # perfect!
# # it stores the past value, 4, skipping the others; at least, that was my desired beahaviour
# 
# cat(crayon::blue("- orphan value test\n"))
# # case of orphan value which was in a bigger cluster
# cl_old = c(1,1,1,1,1,2,2,2,2,2,2,2)
# cl_cur = c(2,2,2,2,2,1,1,1,1,1,1,3)
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# 
# cat(crayon::blue("- complex case of cl_old\n"))
# cl_old = c(2,2,2,2,2,2,2,2,2,2,2,2)
# cl_cur = c(2,2,2,2,2,1,1,1,1,1,1,3)
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# 
# 
# cat(crayon::blue("- appeared clusters test\n"))
# cl_old = c(6,6,6,6,5,5,5,5,3,3,3,3,3)
# cl_cur = c(1,1,1,2,2,2,3,3,3,4,4,4,5)
# mode_correct_clusters(cl_old,cl_cur,verbose=1)
# 
# 
# cat(crayon::blue("- another edge case?\n"))
# cl_old = c(1,2,2,3,2,3,3,3,4,3,3,3,3,3,3,3,1,3,1,1,1,3,5,3,3,6,3,3,2,2,2,2,7,2,2,2,8,8,5,9,3,5,5,2,2,5,2,8,3,2,3,1,3,3,3,8,8,5,2,2,3,3,3,8,1,1,1,2,1,2,1,8,2,1,2,2,8,2,2,1,5,3,1,1,3,5,8,3,2,8,3,10,2,2,2,2,2,2,2,5,10,1,1,2,5)
# cl_cur = c(1,8,8,9,2,3,3,3,8,3,8,3,3,3,3,3,3,3,1,4,1,3,5,3,3,3,3,3,2,8,8,2,3,2,8,6,8,8,5,3,3,5,5,2,2,7,2,8,3,2,3,1,3,1,3,8,8,7,2,2,3,3,8,8,1,1,1,2,4,2,1,8,2,1,2,2,10,2,2,1,7,8,3,3,8,7,8,8,2,8,8,10,2,2,2,2,2,2,2,7,10,9,3,2,7)
# mode_correct_clusters(cl_old,cl_cur,very_verbose = 0,verbose=1)

################################################################################
### graph plot

library(igraph)
library(gridExtra)
# g <- sample_gnp(10, 3/10)
# g_mst <- mst(g)
# 
# print(g)
# plot(g)
# plot(g_mst,edge.arrow.size = 0.5)
#
# plot(make_graph(c(1, 2, 2, 3, 3, 4, 5, 6), directed = FALSE))
# 
# #            lati:  A->B      B->C      C->D     D->As
# gg = make_graph(c("A", "B", "B", "C", "C", "D","D","A"), directed = FALSE)
# plot(gg)
# plot(mst(gg))
# 
# get.edgelist(gg)
# get.edgelist(mst(gg))


# si chiama my_sites sennò dava conflitto con sites dei plot :/
my_sites = data.frame(
	id = unique(df_weekly$IDStations),
	x = unique(df_weekly$Longitude),
	y = unique(df_weekly$Latitude)
)

stations = unique(df_weekly$IDStations)

graph_from_string <- function(x) {
	e <- str2expression(strsplit(x, ",")[[1]])
	do.call(igraph:::graph_from_literal_i, list(e))
}

stations_distance = function(st_1,st_2){
	coords_1 = my_sites[which(my_sites$id == st_1),c("x","y")]
	coords_2 = my_sites[which(my_sites$id == st_2),c("x","y")]
	dist2 = (coords_1$x - coords_2$x)^2 + (coords_1$y - coords_2$y)^2
	return(dist2)
}

cat(crayon::red("- assemble_edges(clusters)\n"))
assemble_edges = function(clusters,need_to_debug=0){
	edges_temp_list = NULL
	edges = data.frame(
		x = c(),
		y = c(),
		xend = c(),
		yend = c(),
		cluster = c()
	)
	# for(cl in 1:max(clusters)){
	for(cl in unique(clusters)){
		stations_here = stations[which(clusters==cl)]
		
		# graph_from_literal( A:B:C:D -- A:B:C:D )
		# crea un grafo con tutte le connessioni
		x = paste(paste(stations_here, collapse = ":"),"--",paste(stations_here, collapse = ":"))
		gg=graph_from_string(x)
		# plot(gg)
		# plot(mst(gg))
		
		# define pesi
		elist = get.edgelist(gg)
		pesi=c()
		if(need_to_debug==1){
			cat("dealing with cluster",cl,"\n")
		}
		if(size(elist)[1]>0){
			for(i in 1:size(elist)[1]){
				pesi = c(pesi,stations_distance(elist[i,1],elist[i,2]))
			}
			elist = get.edgelist(mst(gg,weights = pesi))
		} else{
			elist = get.edgelist(mst(gg))
		}
		
		edges_temp = data.frame(
			x = c(),
			y = c(),
			xend = c(),
			yend = c(),
			cluster = c()
		)
		if(size(elist)[1]>0){
			for (i in 1:size(elist)[1]){
				# c'è un problema con STA- nei nomi delle stazioni, che vengono presi come STA lato ecc
				edges_temp[i,c("x","y")] = my_sites[which(my_sites$id==elist[i,1]),c("x","y")]
				edges_temp[i,c("xend","yend")] = my_sites[which(my_sites$id==elist[i,2]),c("x","y")]
				edges_temp[i,"cluster"] = cl
			}
		} else {
			edges_temp[1,c("x","y")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,c("xend","yend")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,"cluster"] = cl
		}
		edges_temp_list[[cl]] = edges_temp
		edges = rbind(edges,edges_temp)
		# cat("current size",size(edges),"\n")
	}
	edges$cluster = factor(edges$cluster)
	return(edges)
}

cat(crayon::red("- assemble_edges_list(clusters)\n"))
assemble_edges_list = function(clusters,need_to_debug=0){
	edges_temp_list = NULL
	edges = data.frame(
		x = c(),
		y = c(),
		xend = c(),
		yend = c(),
		cluster = c()
	)
	# for(cl in 1:max(clusters)){ # non max perché alcuni cluster potrebbero sparire
	# nel senso che ci possono essere salti, dati dalla cluster mode correction
	# tipo 1 2 3 5 6, manca il cluster 4. Ma ora con unique funziona perché evita quelli mancanti
	for(cl in unique(clusters)){
		stations_here = stations[which(clusters==cl)]
		
		# graph_from_literal( A:B:C:D -- A:B:C:D )
		# crea un grafo con tutte le connessioni
		x = paste(paste(stations_here, collapse = ":"),"--",paste(stations_here, collapse = ":"))
		gg=graph_from_string(x)
		# plot(gg)
		# plot(mst(gg))
		
		# define pesi
		elist = get.edgelist(gg)
		pesi=c()
		if(need_to_debug==1){
			cat("dealing with cluster",cl,"\n")
		}
		if(size(elist)[1]>0){
			for(i in 1:size(elist)[1]){
				pesi = c(pesi,stations_distance(elist[i,1],elist[i,2]))
			}
			elist = get.edgelist(mst(gg,weights = pesi))
		} else{
			elist = get.edgelist(mst(gg))
		}
		
		edges_temp = data.frame(
			x = c(),
			y = c(),
			xend = c(),
			yend = c(),
			cluster = c()
		)
		if(size(elist)[1]>0){
			for (i in 1:size(elist)[1]){
				# c'è un problema con STA- nei nomi delle stazioni, che vengono presi come STA lato ecc
				edges_temp[i,c("x","y")] = my_sites[which(my_sites$id==elist[i,1]),c("x","y")]
				edges_temp[i,c("xend","yend")] = my_sites[which(my_sites$id==elist[i,2]),c("x","y")]
				edges_temp[i,"cluster"] = cl
			}
		} else {
			edges_temp[1,c("x","y")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,c("xend","yend")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,"cluster"] = cl
		}
		edges_temp_list[[cl]] = edges_temp
		edges = rbind(edges,edges_temp)
		# cat("current size",size(edges),"\n")
	}
	edges$cluster = factor(edges$cluster)
	return(edges_temp_list)
}


# test zone
# time = 4
# df_cluster_cut = df_cluster[df_cluster$Time==time,]
# clusters = df_cluster_cut$clusters
# edges = assemble_edges(clusters)

cat(crayon::bold("plotter library required!\n"))




cat(crayon::red("- plot_intensities(df_cluster_cut)\n"))
plot_intensities = function(df_cluster_cut,verbose=F) {
	clust_vals = df_cluster_cut$clusters[1:105]
	ycurrent = y[,paste0("w",time)]
	
	new_cols = colora(105,111,0)
	ord = order(ycurrent,decreasing = T)
	if(verbose==T){
		print(tail(ord))
	}
	
	
	par(mar=c(0,0,0,0),pty="s")
	plot(df_cluster_cut$Longitude[ord],df_cluster_cut$Latitude[ord],col=new_cols,pch=19)
	text(df_cluster_cut$Longitude,df_cluster_cut$Latitude-0.04,labels=1:105,cex=0.6,col="gray")
	
	
	# titolo=paste("Cluster map - time",time)
	# 
	# p = ggplot() +theme_bw()+
	# 	geom_sf(data = altre_regioni,
	# 			# fill = color_empty ,color = color_fill,
	# 			linewidth = 0.1,alpha=0.1, show.legend = FALSE) +
	# 	coord_sf(xlim = range(sites$longitude) + padding, ylim = range(sites$latitude) + padding, expand = FALSE)+
	# 	geom_point(data = df_temp, aes(x = Longitude, y = Latitude,
	# 								   # color = cols[factor(clusters)]), size = 2)+
	# 								   # color = new_cols[clusters_now]
	# 								   # color = lerp_pm10_color(ycurrent,colmin="yellow",colmax = "#990000")
	# 								   color = as.list(new_cols)
	# 								   # fill = "white"
	# 								   ),
	# 			   # size = 2
	# 			   size=lerp_pm10_radius(ycurrent,rmax=2)
	# 			   )+theme(legend.position = "none")+labs(title = titolo)

}




cols_default=colora(12,"div",0)[-2]
xlims = c(-2.5,1.5)
cat(crayon::red("- get_graph_plot(df_cluster_cut)\n"))
get_graph_plot = function(df_cluster_cut,cols=cols_default,
						  titolo=paste("Cluster map - time",time),verbose=0,legenda=1){ # already mode_corrected
	clusters_now = df_cluster_cut$clusters
	edges_list = assemble_edges_list(clusters_now)
	
	p  <-  DefaultPlot(add_bg=FALSE)+
		geom_point(data = df_cluster_cut, aes(x = Longitude, y = Latitude,
											  # color = cols[factor(clusters)]), size = 2)+
											  color = cols[clusters_now]), size = 2)+
		labs(title = titolo)
	
	q = p
	actual_clusters = c()
	for(cl in 1:len(edges_list)){
		edges_to_plot = edges_list[[cl]]
		if(!is.null(edges_to_plot)){
			actual_clusters = c(actual_clusters,cl)
			if(verbose==1){
				cat("procesing cluster",cl,"\n")
				print(actual_clusters)
			}
		q = q + geom_segment(aes(x = x, y = y, xend = xend, yend = yend,
								 color = cols[as.numeric(cluster)]),
							 # color = paste0("cl",cl)),
							 linewidth=1.2,
							 data = edges_to_plot,show.legend = TRUE)+
			guides(color = guide_legend(title = "Clusters"))
		}
	}
	q = q +
		scale_colour_identity(guide="legend",labels=paste0("cl",actual_clusters),
							  breaks=cols[actual_clusters])
	if(legenda==0){
		q = q+theme(legend.position = "none")
	}
	return(q)
}


# needs palette from cold to warm
cat(crayon::red("- color_correct_clusters(df_cluster_cut)\n"))
color_correct_clusters = function(df_cluster_cut,idea=1,verbose=0,max_overall_clusters=8,nint=15) {
	clusters_now = df_cluster_cut$clusters
	palette_heat = 111 # o 77 per ora le migliori
	
	# FUN = mean
	media = c()
	ycurrent = y[,paste0("w",time)]
	cls_labels = unique(clusters_now)
	for (i in 1:length(cls_labels)) {
		media[i] = mean(ycurrent[which(clusters_now==cls_labels[i])])
	}
	res <- data.frame(names=cls_labels, media=media)
	ord = order(res$media,decreasing=T)
	# res[ord,]
	# cat("order =",ord,"\n")
	if(verbose==1){
	for (i in 1:length(cls_labels)) {
		media[i] = mean(ycurrent[which(clusters_now==cls_labels[i])])
		cat(paste0("Media cl",cls_labels[i], " = ",round(media[i],4)
				   ," [pos ",which(ord==cls_labels[i]),"]\n"))
	}
	}
		# clusters_now_cols_ord = clusters_now
	# for (cl in cls_labels){
	# 	clusters_now_cols_ord[which(clusters_now_cols_ord==cl)] = letters[which(ord==cl)]
	# }
	# clusters_now_cols_ord_num = rep(0,length(clusters_now_cols_ord))
	# for (i in 1:length(clusters_now_cols_ord)){
	# 	clusters_now_cols_ord_num[i] = which(letters==clusters_now_cols_ord[i])
	# }
	# cat(clusters_now,"\n")
	# cat(clusters_now_cols_ord_num,"\n")
	if(idea==1){
	########### original idea
	max_overall_clusters = max_overall_clusters
	cols_original = colora(max_overall_clusters,palette_heat,0)

	cols = cols_original
	for (i in 1:length(cls_labels)){
		cols[ord[i]] = cols_original[i] # origial idea
		# 	# left: ord res
		# 	#                right: ordine colori
		# 	cols[6] = cols_original[1]
		# 	cols[2] = cols_original[2]
		# 	cols[3] = cols_original[3]
		# 	cols[5] = cols_original[4]
		# 	cols[4] = cols_original[5]
		# 	cols[1] = cols_original[6]
		# 	cols[7] = cols_original[7]
	}
	}
	if(idea==2){
	############ grid idea
	# a = min(-1,min(media)); b=max(1,max(media));
	a = -1.8; b=1.2;
	nint = nint 
	griglia = seq(a,b,length.out = nint)
	
	cols_original = rev(colora(nint,palette_heat,0)) # if grid idea
	cols = cols_original
	# for (i in 1:length(cls_labels)){
	# for (i in order(media,decreasing = T)){
	for (i in ord){
		target_index = which.min(abs(griglia-media[i]))
		cols[res$names[i]] = cols_original[target_index] # grid idea
		if(griglia[target_index]==1000) {
			# cat(crayon::silver("Tied colors. Suggestion: increment nint\n"))
			stop("Tied colors. Suggestion: increment nint\n")
		}
		# al cluster i va il colore in posizione più vicino al suo valore di media nella griglia
		if(verbose==1){
		cat("cluster",i,"pos griglia/color index",target_index,"\n")
		}
		griglia[target_index:nint] = 1000 # per rimuovere casi di pareggio colore
	}
	}
	return(cols)
}




cat(crayon::red("- get_hist_fill_plot(df_cluster_cut)\n"))
get_hist_fill_plot = function(df_cluster_cut,titolo=paste("Time",time),verbose=0){
	clusters_now = df_cluster_cut$clusters # needs to be already mode corrected if wanted
	# n_clusters = max(clusters_now)
	n_clusters = unique(clusters_now)
	ycurrent = y[,paste0("w",time)]
	
	if(verbose==1){
	cat(crayon::red("Time",time,"\n"))
		# for (cl in 1:n_clusters){
		for (cl in n_clusters){
			cat("Cluster",cl,"- size",length(ycurrent[which(clusters_now==cl)]),
				"- mean",mean(ycurrent[which(clusters_now==cl)]),"\n")
			
		}
	}
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	p = ggplot(df_temp, aes(ycurrent,
							fill = cols[clust_vals]
							# color = cols[clust_vals]
	))+
		
		geom_histogram(alpha=0.3,
					   # fill="white",
					   position="identity")+ # to have histograms
		# geom_density(alpha = 0.3)+ # to have the kernel/density estimation
		
		ggtitle(titolo)+
		# labs(title = paste("Cluster map - time",time))+
		guides(fill = guide_legend(title = "Clusters"))+
		
		# theme_classic()
		theme_bw()+
		# xlim(extrema(ycurrent)+c(-pad,pad))+
		
		# scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
							# breaks=cols[1:max(clust_vals)])+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),
							breaks=cols[n_clusters])+
		# scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
		# breaks=cols[1:max(clust_vals)])
		xlab("log(PM10) values")+
		xlim(xlims)
	
}


cat(crayon::red("- get_hist_color_plot(df_cluster_cut)\n"))
get_hist_color_plot = function(df_cluster_cut,titolo=paste("Time",time),verbose=0){
	clusters_now = df_cluster_cut$clusters # needs to be already mode corrected if wanted
	# n_clusters = max(clusters_now)
	n_clusters = unique(clusters_now)
	ycurrent = y[,paste0("w",time)]
	
	if(verbose==1){
	cat(crayon::red("Time",time,"\n"))
		# for (cl in 1:n_clusters){
		for (cl in n_clusters){
			cat("Cluster",cl,"- size",length(ycurrent[which(clusters_now==cl)]),
				"- mean",mean(ycurrent[which(clusters_now==cl)]),"\n")
			
		}
	}
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	p = ggplot(df_temp, aes(ycurrent,
							# fill = cols[clust_vals]
							color = cols[clust_vals]
	))+
		
		geom_histogram(alpha=0.5,
					   fill="white",
					   position="identity")+ # to have histograms
		# geom_density(alpha = 0.3)+ # to have the kernel/density estimation
		
		ggtitle(titolo)+
		# labs(title = paste("Cluster map - time",time))+
		guides(color = guide_legend(title = "Clusters"))+
		
		# theme_classic()
		theme_bw()+
		# xlim(extrema(ycurrent)+c(-pad,pad))+
		
		# scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
		# breaks=cols[1:max(clust_vals)])
		# scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
							 # breaks=cols[1:max(clust_vals)])+
		scale_color_identity(guide="legend",labels=paste0("cl",n_clusters),
							 breaks=cols[n_clusters])+
		xlab("log(PM10) values")+
		xlim(xlims)
}


cat(crayon::red("- get_hist_continuos_plot(df_cluster_cut)\n"))
get_hist_continuos_plot = function(df_cluster_cut,titolo=paste("Time",time),verbose=1){
	clusters_now = df_cluster_cut$clusters
	# n_clusters = max(clusters_now)
	n_clusters = unique(clusters_now)
	ycurrent = y[,paste0("w",time)]
	cat(crayon::red("Time",time,"\n"))
	
	if(verbose==1){
		# for (cl in 1:n_clusters){
		for (cl in n_clusters){
			cat("Cluster",cl,"- size",length(ycurrent[which(clusters_now==cl)]),
				"- mean",mean(ycurrent[which(clusters_now==cl)]),"\n")
			
		}
	}
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clust_vals=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	# p = ggplot(df_temp, aes(ycurrent, fill = factor(clust_vals))) +
	p = ggplot(df_temp, aes(ycurrent, fill = cols[clust_vals])) +
		
		# scale_fill_manual(values = colora(n_clusters,77,0),name="Cluster")+
		# scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
							# breaks=cols[1:max(clust_vals)])+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),
							breaks=cols[n_clusters])+
		guides(fill = guide_legend(title = "Clusters"))+
		
		geom_density(alpha = 0.3)+
		ggtitle(titolo)+
		theme_bw()+
		xlab("log(PM10) values")+
		ylab("")+
		xlim(xlims)
}


cat(crayon::red("- get_boxplot_plot(df_cluster_cut)\n"))
get_boxplot_plot = function(df_cluster_cut,cols=cols_default,titolo=paste("Time",time),annotate=FALSE){
	clusters_now = df_cluster_cut$clusters # needs to be already mode corrected if wanted
	# n_clusters = max(clusters_now)
	n_clusters = unique(clusters_now)
	ycurrent = y[,paste0("w",time)]
	clsize = table(clusters_now)
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	p = ggplot(df_temp, aes(as.factor(clusters),ycurrent,
							fill = cols[clusters]
							# color = cols[clust_vals]
	))+
		geom_boxplot()+
		# geom_jitter(width=0.2)+
		ggtitle(titolo)+
		labs(title = paste("Cluster map - time",time))+
		guides(fill = guide_legend(title = "Clusters"))+
		
		# theme_classic()
		theme_bw()+
		xlab("clusters")+
		ylab("log(PM10) values")+
		ylim(xlims)+
		# xlim(extrema(ycurrent)+c(-pad,pad))+
		# scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
		# breaks=cols[1:max(clust_vals)])+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),
							breaks=cols[n_clusters])
	if(annotate==TRUE){
		p = p+ annotate("text", x = n_clusters, y = 1.3, label = paste0("size: ",as.vector(clsize)),col="gray")
	}
	return(p)
	# scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
	# breaks=cols[1:max(clust_vals)])
}

cat(crayon::red("- get_boxplot_covariate_plot(df_cluster_cut)\n"))
get_boxplot_covariate_plot = function(df_cluster_cut,cols=cols_default,titolo=paste("Time",time),annotate=FALSE,covariata="AQ_pm10",global_ylim=T){
	clusters_now = df_cluster_cut$clusters # needs to be already mode corrected if wanted
	# n_clusters = max(clusters_now)
	n_clusters = unique(clusters_now)
	
	# ycurrent = y[,paste0("w",time)]
	times = unique(df_wsc$Time)
	# cat("using time=",time,"\n")
	df_time_t = df_wsc[df_wsc$Time == times[time],covariata]
	ycurrent = df_time_t
	
	clsize = table(clusters_now)
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	# print(df_temp)
	
	if (global_ylim == T){
	ylims = c(min(df_wsc[df_wsc$Time<=times[12],covariata]),max(df_wsc[df_wsc$Time<=times[12],covariata]))
	} else {
		ylims = c(min(ycurrent),max(ycurrent))
	}
	
	pad = 2	
	p = ggplot(df_temp, aes(as.factor(clusters),ycurrent[[1]],
							fill = cols[clusters]
							# color = cols[clust_vals]
	))+
		geom_boxplot(show.legend = FALSE)+
		# geom_jitter(width=0.2)+
		ggtitle(titolo)+
		# labs(title = paste("Cluster map - time",time))+
		# labs(title = paste(covariata,"- time",time))+
		labs(title = covariata)+
		guides(fill = guide_legend(title = "Clusters"))+
		
		# theme_classic()
		theme_bw()+
		xlab("")+
		# ylab("log(PM10) values")+
		# ylab(covariata)+
		ylab("")+
		ylim(ylims)+
		# xlim(extrema(ycurrent)+c(-pad,pad))+
		# scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
		# breaks=cols[1:max(clust_vals)])+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),
							breaks=cols[n_clusters])
	if(annotate==TRUE){
		p = p+ annotate("text", x = n_clusters, y = 1.3, label = paste0("size: ",as.vector(clsize)),col="gray")
	}
	return(p)
	# scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
	# breaks=cols[1:max(clust_vals)])
}


library(gridExtra)
cat(crayon::red("- plot_graph_and_hist(df_cluster_cut)\n"))
plot_graph_and_hist = function(df_cluster_cut,cols=cols_default,titolo="",annotate=FALSE,jittera=FALSE){
clusters_now = df_cluster_cut$clusters
# GRAPH #######################
q_graph = get_graph_plot(df_cluster_cut,cols,titolo = paste(titolo,"// Time",time))
# HIST #######################
# by hand as we have to remove the legend here, while the function produces it
# n_clusters = max(clusters_now)
n_clusters = unique(clusters_now)
clsize = table(clusters_now)
ycurrent = y[,paste0("w",time)]
clust_vals = clusters_now[1:105]
df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)


########## HIST
# p = ggplot(df_temp, aes(ycurrent,
# 						fill = cols[clust_vals] # case FILL
# 						# color = cols[clust_vals] # case COLOR
# ))+
# 	geom_histogram(alpha=0.5,
# 				   # fill="white", # case COLOR
# 				   position="identity")+
# 	ggtitle(titolo)+
# 	# guides(color = guide_legend(title = "Clusters"))+
# 	theme_bw()+
# 	theme(legend.position = "none")+
# 	# scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
# 						 # breaks=cols[1:max(clust_vals)])+
# 	scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters), # case FILL
# 	# scale_color_identity(guide="legend",labels=paste0("cl",n_clusters), # case COLOR
# 						 breaks=cols[n_clusters])+
# 	xlab("log(PM10) values")+
# 	xlim(xlims)
# p

############ BOXPLOT
p = ggplot(df_temp, aes(as.factor(clusters),ycurrent,
						fill = cols[clusters]
						# color = cols[clust_vals]
))+
	geom_boxplot(position = "identity")
if(jittera==TRUE){
	p = p+ geom_jitter(width=0.2,size=1,col="#9202af")
}
p = p +
	theme_bw()+
	theme(legend.position = "none")+
	scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),
						breaks=cols[n_clusters])+
	ylab("log(PM10) values")+
	xlab("clusters")+
	ylim(xlims)

	if(annotate==TRUE){
		p = p+ annotate("text", x = n_clusters, y = 1.3, label = paste0("size\n",as.vector(clsize)),col="gray")
	}


p = grid.arrange(q_graph, p, ncol=2,widths=c(1.8,1.2))
# p = arrangeGrob(q_graph, p, ncol=2,widths=c(1.8,1.2))
}



Federica_covariates_plot = function(df_cluster_cut,cols=cols_default,titolo=paste("Time",time),covariates_idx){
	clusters_now = df_cluster_cut$clusters
	
	# GRAPH #######################
	q_graph = get_graph_plot(df_cluster_cut,cols,titolo = titolo)
	n_clusters = unique(clusters_now)
	clsize = table(clusters_now)
	ycurrent = y[,paste0("w",time)]
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	
	############ BOXPLOT PM10
	p_boxplot_pm = ggplot(df_temp, aes(as.factor(clusters),ycurrent,
							fill = cols[clusters]
							# color = cols[clust_vals]
	))+geom_boxplot(position = "identity") +
	theme_bw()+
	theme(legend.position = "none")+
	scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),
						breaks=cols[n_clusters])+
	ylab("log(PM10) values")+
	xlab("clusters")+
	ylim(xlims)
	
	
	############ BOXPLOT covariate
	df_wsc_week = df_wsc[which(df_wsc$week==time),]
	df_temp$Altitude = df_wsc_week$Altitude
	df_temp$WE_wind_speed_100m_max = df_wsc_week$WE_wind_speed_100m_max
	df_temp$LA_lvi = df_wsc_week$LA_lvi
	df_temp$EM_nox_sum = df_wsc_week$EM_nox_sum

										   
	p_boxplot_Altitude = ggplot(df_temp, aes(as.factor(clusters),
											 Altitude,fill = cols[clusters]))+
		geom_boxplot(position = "identity") +
		theme_bw()+
		theme(legend.position = "none")+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),breaks=cols[n_clusters])+
		ylab("Altitude")+
		xlab("clusters")+
		ylim(extrema(df_wsc$Altitude))
	
	p_boxplot_WE_wind_speed_100m_max = ggplot(df_temp, aes(as.factor(clusters),
														   WE_wind_speed_100m_max,fill = cols[clusters]))+
		geom_boxplot(position = "identity") +
		theme_bw()+
		theme(legend.position = "none")+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),breaks=cols[n_clusters])+
		ylab("WE_wind_speed_100m_max")+
		xlab("clusters")+
		ylim(extrema(df_wsc$WE_wind_speed_100m_max))
	
	p_boxplot_LA_lvi = ggplot(df_temp, aes(as.factor(clusters),
										   LA_lvi,fill = cols[clusters]))+
		geom_boxplot(position = "identity") +
		theme_bw()+
		theme(legend.position = "none")+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),breaks=cols[n_clusters])+
		ylab("LA_lvi")+
		xlab("clusters")+
		ylim(extrema(df_wsc$LA_lvi))
	
	p_boxplot_EM_nox_sum = ggplot(df_temp, aes(as.factor(clusters),
											   EM_nox_sum,fill = cols[clusters]))+
		geom_boxplot(position = "identity") +
		theme_bw()+
		theme(legend.position = "none")+
		scale_fill_identity(guide="legend",labels=paste0("cl",n_clusters),breaks=cols[n_clusters])+
		ylab("EM_nox_sum")+
		xlab("clusters")+
		ylim(extrema(df_wsc$EM_nox_sum))
	
	
	library(lubridate)
	# data_from_to = data_agc_lomb[intersect(which(data_agc_lomb$Time>=as.Date("2018-01-01")),
								 	# which(data_agc_lomb$Time<=as.Date("2018-12-31"))),]
	# data_from_to$week = week(data_from_to$Time)
	# data_from_to = data_from_to[which(data_from_to$week==time),]
	# OR this
	data_from_to = df_wsc_week
	
	wind_arrows <- data.frame(
		longitude = sites_plt$Longitude,
		latitude = sites_plt$Latitude,
		direction = cardinal_to_degree(data_from_to$WE_mode_wind_direction_100m),
		Intensity = as.numeric(data_from_to$WE_wind_speed_100m_mean)
	)
	# Calcola le coordinate di fine delle frecce in base alla direzione 
	# l'intensità verra colorata invece di cambiare in lunghezza
	wind_arrows$end_longitude =( wind_arrows$longitude + sin(wind_arrows$direction)/10)
	wind_arrows$end_latitude = ( wind_arrows$latitude + cos(wind_arrows$direction)/10)
	# put a 0 in the NA
	wind_arrows[is.na(wind_arrows)] <- 0
	mappa_wind <- ggplot(data = wind_arrows) +
		# geom_sf(data = shp_map, fill = color_background_map , color = "black", linewidth = 0.5,alpha=0.6)+
		geom_sf(data = altre_regioni, fill = color_empty ,color = color_fill, linewidth = 0.1,alpha=0.1, show.legend = FALSE) +
		geom_sf(data = lombardia_2, fill = color_empty, color = color_comuni_lombardia, linewidth = 0.3,alpha=0.7, show.legend = FALSE) +
		
		coord_sf(xlim = range(sites$longitude) + padding, ylim = range(sites$latitude) + padding, expand = FALSE)+
		# coord_sf(xlim = range(na.omit(wind_arrows$longitude))+padding,
				 # ylim = range(na.omit(wind_arrows$latitude))+padding, expand = FALSE)+
		geom_segment(data = wind_arrows,
					 aes(x = longitude, y = latitude, xend = end_longitude, yend = end_latitude,
					 	color = Intensity),
					 arrow = arrow(angle=10, type = "open", length = unit(0.2, "inches"), ends = "last"),
					 lineend = "round", linewidth = 0.4,alpha=0.9 )+
		theme_bw()+
		theme(panel.grid = element_blank())+
		labs(title="Wind map")+
		ylab("")+
		xlab("")
	
	p = grid.arrange(q_graph,
					 p_boxplot_pm,
					 p_boxplot_Altitude,
					 p_boxplot_WE_wind_speed_100m_max,
					 p_boxplot_LA_lvi,
					 p_boxplot_EM_nox_sum,
					 mappa_wind,
					 ncol=2,
					 # widths=c(1.8,1.2),
					 layout_matrix = rbind(c(1,1,1,1,2,2,2),
					 					   c(1,1,1,1,2,2,2),
					 					   c(1,1,1,1,2,2,2),
					 					   # c(1,1,1,1,2,2),
					 					   c(3,4,5,6,7,7,7),
					 					   c(3,4,5,6,7,7,7))
	)
	return(p)
}


easy_plot = function(clusters_input,nintput=30){
	clusters_input = clusters_input[1:105]
	sites = data.frame(
		longitude = unique(df_weekly$Longitude), 
		latitude = unique(df_weekly$Latitude))
	std_sites = data.frame(
		longitude = unique(df_wsc$Longitude), 
		latitude = unique(df_wsc$Latitude))
	stations = unique(df_wsc$IDStations)
	y=data.frame()
	for(st in stations){
		y_we_pm10=cbind(as.data.frame(st),t(df_wsc[which(df_wsc$IDStations==st),"AQ_pm10"]))
		y=rbind(y,y_we_pm10)
	}
	rownames(y) = NULL
	colnames(y)<- c("id",paste0("w", 1:53))
	
	df_temp = data.frame(
		Longitude = sites$longitude,
		Latitude = sites$latitude,
		clusters = clusters_input
	)
	cols = color_correct_clusters(df_temp,idea=2,verbose=0,nint=nintput)
	p = plot_graph_and_hist(df_temp,cols,jittera = T)
}
