# Update the Beta Paramter and Recalculate the relevant internal matrices
#
#
library(universalmotif)
library(memes)
options(meme_bin = "/share/pkg.7/meme/5.3.3/install/bin/") #Set the meme bin directory

change_beta_update_ppm <- function(corec_motif,beta,name="motif") {

    zscore = seed_zscore(corec_motif@zscore_motif) #Uses Rebekah's seed_zscore function (motif_scores.R)

    if(beta == "linear"){

           if(zscore > 6) { corec_motif@beta = 1}  #Constant above cutoff
      else if(zscore < 0) { corec_motif@beta = 4} #Constant below zero (will probably disregard these)
      else                { corec_motif@beta = 4 - 3/6*zscore} #Linear from zero to 6

     }else{
      corec_motif@beta = beta
    }

    # Transform the z-scores using the beta parameter

    corec_motif@ppm  <-
        # Multiply each z-score by beta and then take the exponential
        exp(corec_motif@beta * corec_motif@zscore_motif) %>%

        # Normalize by dividing each value in each column by the column-wise sum
        # dplyr::mutate_all(~ . / sum(.))  -- also works
        # dplyr::mutate_all(.fun = ~ . / sum(.))) -- also works
        dplyr::mutate_all(list(~ . / sum(.))) %>%

        # Convert to a numeric matrix
        as.matrix() %>%

        # Convert to a universalmotif object
        universalmotif::create_motif(name = name)

    # Return the corec motif
    return(corec_motif)
}


plot_universal_motif <-
    function(
        reference_motif,
        output_logo_type = c("PPM", "PWM", "ICM"),
        pval = 1
    ) {

        # Make sure the selected reference motif logo type is a valid option
        output_logo_type <-
            match.arg(output_logo_type)

        # Get the correct form of the reference motif as a numeric matrix
        motif_matrix <-
            # Get the correct form of the motif
            switch(
                output_logo_type,
                "PPM" = reference_motif@motif,
                "PWM" =
                    universalmotif::convert_type(
                        reference_motif@motif,
                        "PWM"
                    )@motif,
                "ICM" =
                    universalmotif::convert_type(
                        reference_motif,
                        "ICM"
                    )@motif
            ) %>%

            # Convert to a numeric matrix
            as.matrix()

        # Figure out the width of the reference motif
        motif_width <-
            ncol(motif_matrix)

        # Begin the plot of the reference motif that matches the corecmotif
        motif_plot <-
            ggseqlogo::ggseqlogo(
               motif_matrix,
                method = "custom",
                seq_type = "dna"
            )

        # Get the range of the y axis
        motif_plot_yrange <-
            ggplot2::ggplot_build(
                motif_plot
            )$layout$panel_params[[1]]$y.range

        # Format the width and height of the axes without printing any messages
        suppressMessages(
            motif_plot <-
                motif_plot +

                # Increase the plot width by decreasing the padding on the edges
                # ggseqlogo already makes an x scale, so this prints a message
                ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +

                # Increase the plot height and include 1 decimal place in labels
                ggplot2::scale_y_continuous(
                    expand = c(0.01, 0.01),
                    labels = function(x) {sprintf("%.1f", x)}
                )
        )

        # Figure out what the y axis label should be based on the motif type
        y_label <-
            switch(
                output_logo_type,
                "PPM" = "probability",
                "PWM" = "weight",
                "ICM" = "bits"
            )

        # Add remaining plot elements
        motif_plot <-
            motif_plot +

            # Add a green outline to the logo
            ggplot2::annotate(
                'rect',
                xmin = 0.5,
                xmax = motif_width + 0.5,
                ymin = motif_plot_yrange[1],
                ymax = motif_plot_yrange[2],
                color = "#added1",
                fill = NA,
                size = 1.5
            ) +

            # Add a title with the name of the reference motif and the logo type
            ggplot2::ggtitle(
                paste(
                    "Motif",
                    paste(
                        reference_motif@name,
                        reference_motif@altname,
                        sep = "_"
                    ),
                    paste(
                        "Match p-value:",
                        signif(pval, 3)
                    ),
                    output_logo_type,
                    sep="\n"
                )
            ) +

            # Add the y axis label
            ggplot2::ylab(y_label) +

            # Set the formatting for the axis text and the plot title
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = 10, face = "plain"),
                axis.title.y = ggplot2::element_text(size = 10, face = "plain"),
                plot.title = ggplot2::element_text(size = 10)
            )


        # Return the final plot
        return(motif_plot)
    }

plot_ppm_matrix <-
    function(
        motif,
        output_logo_type = c("PPM", "PWM", "ICM")
    ) {

        # Make sure the selected reference motif logo type is a valid option
        output_logo_type <-
            match.arg(output_logo_type)

        # Get the correct form of the reference motif as a numeric matrix
        motif_matrix <-
            # Get the correct form of the motif
            switch(
                output_logo_type,
                "PPM" = motif,
                "PWM" =
                    universalmotif::convert_type(
                        motif,
                        "PWM"
                    )@motif,
                "ICM" =
                    suppressMessages(universalmotif::convert_type(
                        motif,
                        "ICM"
                    )@motif)
            ) %>%

            # Convert to a numeric matrix
            as.matrix()

        # Figure out the width of the reference motif
        motif_width <-
            ncol(motif_matrix)

        # Begin the plot of the reference motif that matches the corecmotif
        motif_plot <-
            ggseqlogo::ggseqlogo(
                motif_matrix,
                method = "custom",
                seq_type = "dna"
            )

        # Get the range of the y axis
        motif_plot_yrange <-
            ggplot2::ggplot_build(
                motif_plot
            )$layout$panel_params[[1]]$y.range

        # Format the width and height of the axes without printing any messages
        suppressMessages(
            motif_plot <-
                motif_plot +

                # Increase the plot width by decreasing the padding on the edges
                # ggseqlogo already makes an x scale, so this prints a message
                ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +

                # Increase the plot height and include 1 decimal place in labels
                ggplot2::scale_y_continuous(
                    expand = c(0.01, 0.01),
                    labels = function(x) {sprintf("%.1f", x)}
                )
        )

        # Figure out what the y axis label should be based on the motif type
        y_label <-
            switch(
                output_logo_type,
                "PPM" = "probability",
                "PWM" = "weight",
                "ICM" = "bits"
            )

        # Add remaining plot elements
        motif_plot <-
            motif_plot +

            # Add a green outline to the logo
            ggplot2::annotate(
                'rect',
                xmin = 0.5,
                xmax = motif_width + 0.5,
                ymin = motif_plot_yrange[1],
                ymax = motif_plot_yrange[2],
                color = "#added1",
                fill = NA,
                size = 1.5
            ) +

            # Add a title with the name of the reference motif and the logo type
            ggplot2::ggtitle(
                paste(
                    "Motif",
                    paste(
                        "motif",
                        sep = "_"
                    ),
                    output_logo_type,
                    sep="\n"
                )
            ) +

            # Add the y axis label
            ggplot2::ylab(y_label) +

            # Set the formatting for the axis text and the plot title
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = 10, face = "plain"),
                axis.title.y = ggplot2::element_text(size = 10, face = "plain"),
                plot.title = ggplot2::element_text(size = 10)
            )


        # Return the final plot
        return(motif_plot)
    }


#Convert Reference Motif PPM to ICM for comparison
convert_refs_to_icm <- function(ref_motifs)
{
            ref_icms =
                ref_motifs %>%
                purrr::map(convert_ref_to_icm)

            return(ref_icms)
}
convert_ref_to_icm <- function (ref_motif){

    icm = ref_motif %>%
            universalmotif::convert_type("ICM")

    ref_motif@motif = icm@motif
    return(ref_motif)
}

convert_universal_motif_ppm_to_icm <- function (motif){

  icm = motif %>%
    universalmotif::convert_type("ICM")

  motif@motif = icm@motif
  return(motif)
}

#Conver CoRec Motif PPM to ICM for comparison
#Note - you can reset to PPM by reseting the Beta Paramter (change_beta_update_ppm)
convert_corec_motifs_to_icm <- function(corec_motifs)
{

  #Return a list of CoRec motifs with ICM as the PPM
  corec_icms =
    corec_motifs %>%
    purrr::map(convert_corec_motif_to_icm)

  return(corec_icms)
}
convert_corec_motif_to_icm <- function (corec_motif){

  #Get ICM from PPM
  icm = corec_motif@ppm %>%
    universalmotif::convert_type("ICM")

  #Pass ICM matrix to Corec_Motif@ppm matrix
  corec_motif@ppm@motif = icm@motif
  return(corec_motif)
}

## IDENTIFY_COREC_MOTIF_MATCH_TWS
#
# ml = corec motif
# rm=  reference motifs (list of universal motifs)
#
# NOTE: If you add a new method - need to change code in 2 places (ARGS, and SWITCH)

identify_corec_motif_match_TomTom <-
  function(
    corec_motif,
    database,
    method=c("ed","allr","pearson","sandelin","kullback","blic1","blic5","llr1","llr5"),
    beta = "NA",
    min_overlap=6
  ){

    #TEMP HARD CODE
#    num=33
#    corec_motif = corec_motifs[[num]]
#    database = db
#    min_overlap = 6
#    method = "ed"
#    beta = "linear"
    ###

print(paste0(corec_motif@seed_name))

    # Make sure the select method is allowed
    method <- match.arg(method)

    #Change Beta value if passed
    if(!is.na(beta)) {
        corec_motif = change_beta_update_ppm(corec_motif,beta=beta)
    }

    matches <- runTomTom(corec_motif@ppm,database=database,dist=method, min_overlap=min_overlap) #Run TomTom

    matches_df <-
      matches$tomtom %>%
      as.data.frame()

    if(is.na(matches$best_match_motif)){

      corec_motif@motif_match_score_type <- method
      corec_motif@motif_match_score = 100
      corec_motif@motif_match_pvalue = 1

    }else
    {

      # Update the motif matching slots of the original corecmotif object
      corec_motif@motif_match   = matches_df$match_motif[[1]]
      corec_motif@motif_match_score_type = method
      corec_motif@motif_match_score = matches_df$match_qval[1]
      corec_motif@motif_match_pvalue = matches_df$match_pval[1] #Get top match pvalue
    }#else

    # Return the updated corecmotif
    return(corec_motif)

  } #end function

##
# COMPARE_TWO_UNIVERSAL_MOTIFS
#
# ARGS:  m1, m2 are EACH universalmotif OBJECTS  (i.e., not LISTS)
#        len_norm = TRUE|FALSE
#        ovlerap = NA -> overlap is length of shortest motif
#                  alternatively you set the overlap min (e.g., 8 bp)
#
## NOTE: If you add a new method - need to change code in 3 places in this function (ARGS, and 2 SWITCHES)
compare_two_universal_motifs <-
    function(
        m1,m2,
        method=c(
          "EUCL",
          "PCC",
          "EUCL_ICM",
          "PCC_ICM"),
        len_norm=TRUE,
        overlap="NA"
    ){

## TEMP HARD CODE ARGS
#    num=1
#    m1=test_ref[[4]]
#    m2=corec_motifs[[num]]@ppm
#    plot_ppm_matrix(m1@motif,"PPM")
#    plot_ppm_matrix(m2@motif,"PPM")

#    method="PCC"
#    method="PCC"
#    len_norm = TRUE
#    overlap="NA"

### EDN TEMP

    # Make sure the selected motif logo type is a valid option
    method <- match.arg(method)

    # If Method uses ICM - convert m1 and m2 to ICMs
    if(grepl("ICM",method)){
      m1 = convert_universal_motif_ppm_to_icm(m1)
      m2 = convert_universal_motif_ppm_to_icm(m2)
    }

    #Set LONG and SHORT matrices
    if(ncol(m1@motif) < ncol(m2@motif)) {
        ml = as.matrix(m2@motif) #ml - Matrix Long
        ms_for = as.matrix(m1@motif) #ms - Matrix short

        #Get RC motif
        ms_rc =  as.matrix(
                    return_motif_rc(m1@motif)
                    )
    }else{
        ml = as.matrix(m1@motif)
        ms_for = as.matrix(m2@motif)

        #Get RC motif
        ms_rc =  as.matrix(
                    return_motif_rc(m2@motif)
                    )
    }

    #Make sure shortest motif is not shorter than overlap
    ns = ncol(ms_for)
    nl = ncol(ml)

    # If overlap not set - reset to min motif length
    if(overlap == "NA") { overlap = ns}

    ### TEMP PLOT
#    pl = list()
#    pl[[1]] = plot_ppm_matrix(ml,"PPM")
#    pl[[2]] = plot_ppm_matrix(ms_for,"PPM")
#    pl[[3]] = plot_ppm_matrix(ms_rc,"PPM")
#    cowplot::plot_grid(plotlist = pl,nrow=2)

    #CHECK overlap versus motif lengths
    if( ns < overlap ){
      ovelap = ns
      print(paste("Trying to compare motifs:: one is shorter (N=",ns,") than overlap of",overlap," re-setting overlap to min motif length"))
#        stop(paste("Trying to compare motifs:: one is shorter (N=",ns,") than overlap of",overlap))
    }

    ### RUBRIC for determining the COLUMN SUBSETS to cmopare while 'sliding'
    #
    # ml =     123456789            (Motif LONG)
    # ms =  1234567                 (Motif SHORT, OVERLAP =4)
    #
    # OVERLAP minimum determines the OFFSET/OVERHANGS on either end
    #       -- +++++++++
    # reg = 210123456789   (Basically the position w.r.t. ml motif)
    # ml =     123456789
    # ms =  1234567        (reg= -2, ms_START=4)
    # ms =   1234567       (reg= -1, ms_START=3)
    # ms =    1234567      (reg=  0, ms_START=2)

    # ms =     1234567
    # ms =         1234567
    #
    #     RUBRIC for determining columns to compare::
    #
    #     reg_START = -1*(ns-overlap)+1  (e.g., -1*(7-4-1) = -2)
    #     reg_STOP  = nl-overlap+1       (e.g., 9-4+1 = 6)
    #
    #     ml_START = reg       if ml_START <= 0, ml_START=1
    #     ml_STOP  = reg+ns-1  if ml_STOP > nl, ml_STOP=nl
    #
    #     ms_START =  1       if reg <= 0, ms_START= 2-reg  (e.g., 2--2=4, 2--1=3, 2--0=2,)
    #     ms_STOP  =  ns      if reg+ns > nl, ms_STOP= 1+ nl-reg

    reg_START = -1*(ns-overlap)+1
    reg_STOP  = nl-overlap+1

    scores = vector()
    scores_or = vector()
    scores_reg = vector()
    data = data.frame()

    orientations = c("for","rc") #FORWARD and REVERSE COMPLMENT

    #LOOP through 2 orientations and ALL registers comparing sub matrices
     for( or in orientations){

          if( or == "for"){ ms = ms_for}
          else{ ms = ms_rc} #If comparing RC - use RC matrix

       for(reg in reg_START:reg_STOP){
#       for(reg in reg_START:(reg_START+10)){
        #Set Column Parameters
        ml_START = reg
            if(reg <= 0){ ml_START = 1}
        ml_STOP = reg+ns-1
            if(ml_STOP > nl) { ml_STOP=nl}
        ms_START = 1
            if(reg <=0){ms_START = 2-reg}
        ms_STOP = ns
            if(reg+ns-nl > 0) { ms_STOP = 1+nl-reg}

        #Get subset matrices
        sub_ml = ml[,ml_START:ml_STOP]
        sub_ms = ms[,ms_START:ms_STOP]

        #Calculate Scores between aligned sub-matrices
        reg_score <-
            # Get the correct form of the motif
            switch(
                method,
                "EUCL"     = score_matrices_EUCL(sub_ml,sub_ms,norm=len_norm),
                "EUCL_ICM" = score_matrices_EUCL(sub_ml,sub_ms,norm=len_norm),
                "PCC"      = score_matrices_PCC(sub_ml,sub_ms,norm=len_norm),
                "PCC_ICM"  = score_matrices_PCC(sub_ml,sub_ms,norm=len_norm)
            )

        scores = scores %>% append(reg_score)
        scores_or = scores_or %>% append(or)
        scores_reg = scores_reg %>% append(reg)

        }#for reg
     }# for orientations

    data = cbind(
              scores,
              scores_or,
              scores_reg ) %>%
          as.data.frame()

    data$scores = as.numeric(data$scores) #For some reason becomes string
    data$scores_reg = as.numeric(data$scores_reg)

    #Get the best scoring register/orientation
    best_score <-
      # Get the correct form of the motif
      switch(
        method,
        "EUCL"     = min(data$scores),
        "EUCL_ICM" = min(data$scores),
        "PCC"      = max(data$scores),
        "PCC_ICM"  = max(data$scores)
      )

    ## TEMP -- Plot comparison details
    plot_details=FALSE
    if(plot_details){
    #data[which(data$scores_or == "for"),]
    data_for = data %>%
          filter(grepl("for",scores_or))
    data_rc = data %>%
      filter(grepl("rc",scores_or))

   score_plot <- data_for %>%
        ggplot(aes(x=scores_reg,y=scores))+
        scale_y_continuous(name="scores",
                           limits=c(min(data$scores),max(data$scores))
                           ) +
        scale_x_continuous(name="scores_reg",
                           n.breaks = max(data$scores_reg) - min(data$scores_reg)
                           ) +
        geom_point()+
        geom_point(data = data_rc, aes(x=scores_reg,y=scores),col="red")
    score_plot

    #Plot motifs at different 'REG'
    reg = 14
    l = return_sub_matrices(reg,nl,ns,ml,ms_for)
    sub_ml = l[[1]]
    sub_ms = l[[2]]
    pl = list()
    pl[[1]] = plot_ppm_matrix(sub_ml,"ICM")
    pl[[2]] = plot_ppm_matrix(sub_ms,"ICM")
    cowplot::plot_grid(plotlist = pl,nrow=2)
    }
    ## END TEMP PLOT ROUTINE

    return(best_score)

} # end function


## IDENTIFY_COREC_MOTIF_MATCH_TWS
#
# ml = corec motif
# rm=  reference motifs (list of universal motifs)
#
# NOTE: If you add a new method - need to change code in 2 places (ARGS, and SWITCH)

identify_corec_motif_match_TWS <-
  function(
    corec_motif,
    ref_motifs,
    method=c("EUCL","PCC","EUCL_ICM","PCC_ICM"),
    len_norm=TRUE,
    overlap="NA"
  ){

#    print(paste0(" matching motif=",corec_motif@seed_name))
#    corec_motif = corec_motifs[[1]]
#    ref_motifs = reference_motifs
#    ref_motifs = test_ref
#    method="PCC_ICM"
#    len_norm=TRUE
#    overlap=10

    #Match all Ref motifs against this CoRec Motif (use PPM)
    match_scores =
        ref_motifs %>%
        purrr::map(~compare_two_universal_motifs(.x,corec_motif@ppm,method,len_norm,overlap)) %>%
        unlist()

    #Identify the best matching motif
    best_match <-
      # Get the correct form of the motif
      switch(
        method,
        "EUCL" =     which.min(match_scores),
        "EUCL_ICM" = which.min(match_scores),
        "PCC" =      which.max(match_scores),
        "PCC_ICM" =  which.max(match_scores)
      )

    # Update the motif matching slots of the original corecmotif object
    corec_motif@motif_match <-
      ref_motifs[[best_match]]

    corec_motif@motif_match_score_type <-
      method

    corec_motif@motif_match_score <-
      match_scores[best_match]

    corec_motif@motif_match_pvalue <-
      1

    # Return the updated corecmotif
    return(corec_motif)

  } #end function


### Identify CoRec entry that matches some SNP_id_allele


# Extract sub matrcies based on the REG counter (used in compare_two_universalmotifs)
return_sub_matrices <- function(reg,nl,ns,ml,ms){

  ml_START = reg
    if(reg <= 0){ ml_START = 1}
  ml_STOP = reg+ns-1
    if(ml_STOP > nl) { ml_STOP=nl}
  ms_START = 1
    if(reg <=0){ms_START = 2-reg}
  ms_STOP = ns
    if(reg+ns-nl > 0) { ms_STOP = 1+nl-reg}

  return( list(ml[,ml_START:ml_STOP],ms[,ms_START:ms_STOP]))
}

# Return Motif RC
#   Needed this as universalmotif::motif_rc re-normalizes column weights to ppm
return_motif_rc <-function(m){
    tm =  m[nrow(m):1,ncol(m):1]
    rownames(tm) = c("A","C","G","T")
    return(tm)
}


# EUCLIDEAN DISTANCE
# Function taken from::
# https://bioconductor.org/packages/release/bioc/vignettes/universalmotif/inst/doc/MotifComparisonAndPvalues.pdf
# c1,c2 = matrices
# Allows for length normalization
score_matrices_EUCL <- function(c1,c2,norm=FALSE) {
        s = sum((c1-c2)^2)
        if(norm) { s = s/(ncol(c1))}
        return(s)
}

# PEARSON CORRELATION COEFFICIENT
# Function taken from::
# https://bioconductor.org/packages/release/bioc/vignettes/universalmotif/inst/doc/MotifComparisonAndPvalues.pdf
# c1,c2 = matrices
# Allows for length normalization
score_matrices_PCC <- function(c1,c2, norm=FALSE) {

  n <- length(c1)
  top <- n * sum(c1 * c2) - sum(c1) * sum(c2)
  bot <- sqrt( ( n * sum(c1^2) - sum(c1)^2 ) * ( n * sum(c2^2) - sum(c2)^2 ) )
  s = top/bot
  if(norm) { s = s/(ncol(c1))}
  return(s)
}

#





