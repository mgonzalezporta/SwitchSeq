## functions

readTranscriptExpressionSet=function(gId=gId, infile=infile, cond1=cond1, cond2=cond2) {
  # header
  command=paste("head -n1", infile, sep=" ")
  header=read.table(pipe(command), as.is=T)

  # conditions
  cond1=.format_cond(cond1)
  cond2=.format_cond(cond2)

  # transcript RPKMs
  command=paste("grep", gId, infile, sep=" ")
  data=read.table(pipe(command), col.names=header,
    colClasses=c("NULL", "character", rep("numeric", length(header)-2)))
  rpkms=as.matrix(data[,c(cond1, cond2)])
  rownames(rpkms)=data[,1]

  # recalculate conditions
  cond1=cond1-1
  cond2=seq(from=length(cond1)+1, to=dim(rpkms)[2], by=1)
  conditions=list(cond1=cond1, cond2=cond2)

  # create TranscriptExpressionSet object
  tes=new("TranscriptExpressionSet",
    id=gId,
    conditions=conditions,
    rpkms=rpkms
  )

  return(tes)
}

annotateTranscriptExpressionSet=function() {
  # include tBiotype info
}

.format_cond=function(cond) {
  x=as.numeric(strsplit(cond, "-")[[1]])
  cond=(x[1]-1):(x[2]-1)
  return(cond)
}

###########

## classes
setClass("TranscriptExpressionSet",
  representation = representation(
    id="character",
    rpkms="matrix",
    conditions="list",
    gexp="numeric",
    dominance="numeric",
    .relexp="matrix",
    .scaledexp="matrix"

    biotypes="data.frame",
    .cols="list"
  )
)

## generics
setGeneric("id", function(object) standardGeneric("id"))
setGeneric("id<-", function(object,value) standardGeneric("id<-"))
setGeneric(".calculate_dominance", function(object) standardGeneric(".calculate_dominance"))
setGeneric(".calculate_relexp", function(object) standardGeneric(".calculate_relexp"))
setGeneric(".ratio_second", function(object) standardGeneric(".ratio_second"))
setGeneric(".calculate_scaledexp", function(object) standardGeneric(".calculate_scaledexp"))

## methods
setMethod("initialize", "TranscriptExpressionSet",
  function(.Object, id=id, rpkms=rpkms, conditions=conditions){
    .Object@id=id
    .Object@rpkms=rpkms
    .Object@conditions=conditions
    .Object@gexp=colSums(.Object@rpkms)
    .Object@dominance=.calculate_dominance(.Object)
    .Object@.relexp=.calculate_relexp(.Object)
    .Object@.scaledexp=.calculate_scaledexp(.Object)
    return(.Object)
  })

setMethod("show",
  "TranscriptExpressionSet",
  function(object) {
    cat("Object of class",class(object),"\n")
    cat("  Id:",id(object),"\n")
    cat("  RPKMs: data frame with", dim(object@rpkms)[1], "transcripts and", dim(object@rpkms)[2]-1, "samples\n")
    cat("    transcripts:", rownames(object@rpkms), "\n")
    cat("    samples:", colnames(object@rpkms), "\n")
    cat("  Conditions:\n")
    cat("    condition 1:", object@conditions$cond1, "\n")
    cat("    condition 2:", object@conditions$cond2, "\n")
    cat("  Gene expression:", head(object@gexp), "(...)\n")
    cat("  Dominance:", head(object@dominance), "(...)\n")
  })

setMethod("id", "TranscriptExpressionSet", function(object) object@id)
setReplaceMethod("id",
                 signature(object="TranscriptExpressionSet",
                           value="character"),
                 function(object, value) {
                   object@id=value
                   return(object)
                 })

setMethod(".calculate_dominance", "TranscriptExpressionSet",
    function(object) {
          dominance=apply(object@rpkms, 2, .ratio_second)
          return(dominance)
    })

setMethod(".calculate_relexp", "TranscriptExpressionSet",
    function(object) {
      relexp=t(t(object@rpkms)/object@gexp)
      return(relexp)
    })

setMethod(".ratio_second", "numeric", 
  function(object) {
    x=object
    x=x[order(x, decreasing=T)]
    r=x[2]/x[1]

    return(r)   
  })

setMethod(".calculate_scaledexp", "TranscriptExpressionSet",
    function(object) {
      x=t( t(object@.relexp) * abs( sqrt( object@gexp/max(object@gexp ))))
      x[is.na(x)]=0
      r=x/max(x)

      return(r)
    })