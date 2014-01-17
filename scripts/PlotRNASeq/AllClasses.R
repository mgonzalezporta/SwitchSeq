setClass("TranscriptExpressionSet",
  representation = representation(
    id="character",
    rpkms="matrix",
    conditions="list",
    gexp="numeric",
    dominance="numeric",
    .relexp="matrix",
    .scaledexp="matrix",
    biotypes="data.frame",
    .cols="list",
    significant_events="character"
  )
)
