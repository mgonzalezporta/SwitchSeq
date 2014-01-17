# accessors
setGeneric("id", 
	function(object) standardGeneric("id"))

setGeneric("id<-", 
	function(object,value) standardGeneric("id<-"))

setGeneric("conditions", 
	function(object) standardGeneric("conditions"))

setGeneric("rpkms", 
	function(object) standardGeneric("rpkms"))

setGeneric("gexp", 
	function(object) standardGeneric("gexp"))

setGeneric("dominance", 
	function(object) standardGeneric("dominance"))

setGeneric("biotypes", 
	function(object) standardGeneric("biotypes"))

setGeneric("significant_events",
        function(object) standardGeneric("significant_events"))

## plots
setGeneric("plotStars", 
	function(tes, outfile) standardGeneric("plotStars"))

setGeneric("plotDistr", 
	function(tes, outfile) standardGeneric("plotDistr"))

## internal
setGeneric(".calculate_dominance", 
	function(object) standardGeneric(".calculate_dominance"))

setGeneric(".calculate_relexp", 
	function(object) standardGeneric(".calculate_relexp"))

setGeneric(".ratio_second", 
	function(object) standardGeneric(".ratio_second"))

setGeneric(".calculate_scaledexp", 
	function(object) standardGeneric(".calculate_scaledexp"))

setGeneric(".get_transcript_colors", 
	function(object) standardGeneric(".get_transcript_colors"))

setGeneric(".distr_summary", 
	function(object) standardGeneric(".distr_summary"))

setGeneric(".get_legend_title_starplot", 
	function(object) standardGeneric(".get_legend_title_starplot"))

setGeneric(".plot_boxplots", 
	function(tes, outfile) standardGeneric(".plot_boxplots"))

setGeneric(".subplot_boxplots", 
	function(ldata=data, xlab=xlab, ylab=ylab, type=type, mar=mar, col=col, plot_count=plot_count) standardGeneric(".subplot_boxplots"))

setGeneric(".plot_segments", 
	function(tes, outfile) standardGeneric(".plot_segments"))

setGeneric(".subplot_segments", 
	function(data1=data1, data2=data2, xlab=xlab, ylab=ylab, type=type, mar=mar, col=col, plot_count=plot_count) standardGeneric(".subplot_segments"))

setGeneric(".annotate_data", 
	function(data=data, significant_events=significant_events, biotypes=biotypes) standardGeneric(".annotate_data"))

setGeneric(".expand_xlab",
        function(object) standardGeneric(".expand_xlab"))
