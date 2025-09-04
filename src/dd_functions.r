library(ggplot2)
library(tidyverse)
library(ggpubr)
library(modelr)
library(stringr)

# Returns a darker shade of a given color
# Example: darker.col("blue", 50)

darker.col = function(color, how.much = 30){
  colorRampPalette(c(color, "black"))(100)[how.much]
}

make_clr_similarity <- function(mat, id_col='genepid', method='cosine', drop_columns=c()) {
    # Select the appropriate method for similarity calculation
    if (method == 'Jaccard') {
        # Convert the matrix to binary (presence/absence) for Jaccard similarity
        mat_binary <- mat %>% select(-drop_columns) %>% column_to_rownames(id_col) %>% as.matrix()
        mat_binary[mat_binary > 0] <- 1
        
        # Calculate Jaccard similarity
        sim.mat <- as.matrix(vegdist(mat_binary, method = 'jaccard', binary = TRUE))
    } else if(method == 'bray'){
        sim.mat <- as.matrix(vegdist(mat %>% select(-drop_columns) %>% column_to_rownames(id_col), method))
    }else {
        sim.mat <- as.matrix(proxy::simil(mat %>% select(-drop_columns) %>% column_to_rownames(id_col), method))
    }
    
    # Convert the similarity matrix to a tidy format
    sim.df <- as_tibble(sim.mat, rownames=NA) %>%
        rownames_to_column(id_col) %>%
        pivot_longer(-c(id_col), values_to = paste0('sim.', method), names_to = paste0(id_col, '_2')) 
    
    return(sim.df)
}

# Create a count matrix of genes by geographic location
make_count_matrix <- function(dat, geoCol='city', geneCol='gene'){
    mat <- dat %>%
        select(all_of(geneCol), all_of(geoCol)) %>%
        table %>%
        as.data.frame.matrix
    return(mat)
}

# Filter a count matrix by minimum gene copies and minimum ARGs per site
filter_count_matrix <- function(mat, min_copies_per_arg=1, min_args_per_site=100) {
    return(mat[rowSums(mat)>min_copies_per_arg, colSums(mat)>min_args_per_site])
}

# Create and filter a count matrix from input data
make_mats <- function(dat, geoCol='city', geneCol='gene', min_copies_per_arg=1, min_args_per_site=100, verbose=F){
    
    # Create unfiltered count matrix
    unfiltered.mat <- make_count_matrix(
        dat = dat, 
        geoCol = geoCol,
        geneCol= geneCol
    )

    # Apply filtering criteria
    filtered.mat <- filter_count_matrix(
        mat= unfiltered.mat, 
        min_copies_per_arg=min_copies_per_arg, 
        min_args_per_site=min_args_per_site
    )

    # Optionally print dimensions
    if(verbose){
        print(paste("Dimension of input data:", nrow(dat), ncol(dat)))
        print(paste("Dimension of unfiltered count matrix:", nrow(unfiltered.mat), ncol(unfiltered.mat)))
        print(paste("Dimension of filtered count matrix:", nrow(filtered.mat), ncol(filtered.mat)))
        
    }

    # Return both matrices
    return(list("mat" = unfiltered.mat, "filtered" = filtered.mat))
}


create_res_dd <- function(mat){
    ddMat <- mat %>%
        t %>%
        decostand(method = "hellinger") %>%
        data.frame(check.names = F) %>%
        as.matrix %>% 
        vegdist() %>% 
        as.matrix() %>% 
        as_tibble(rownames = "area") %>% 
        pivot_longer(-area) %>%
        filter(area != name) %>%
        arrange(value) %>%
        mutate(
            area1 = case_when(
                area < name ~ area,
                area > name ~ name
            ),
            area2 = case_when(
                area < name ~ name,
                area > name ~ area
            )
        ) %>%
        select(-name, -area) %>%
        unique()

    return(ddMat)
}

# Create a pairwise geographic distance table from metadata
create_physical_dd <- function(metadata, geoCol){

    # Clean and prepare sampling point coordinates
    sampling.points <- metadata %>%
        mutate(
            city = case_when(
                city == "Santiago" & country == "Spain" ~ "Santiago_ESP",
                TRUE ~ city)
        ) %>%
        select(all_of(geoCol), lon, lat) %>%
        unique() %>% na.omit %>%
        group_by_(geoCol) %>%
        summarise(lon = median(lon, na.rm = T), lat = median(lat, na.rm = T)) %>%
        st_as_sf(coords = c("lon", "lat"), crs = "WGS84")

    # Compute pairwise geographic distances
    sampling.points.dist <- st_distance(sampling.points, sampling.points) %>%
        data.frame

    # Assign row and column names based on geographic identifiers
    colnames(sampling.points.dist) <- rownames(sampling.points.dist) <- sampling.points[[geoCol]]

    # Reshape distance matrix into long format and clean up
    sampling.points.dist <- sampling.points.dist %>% 
        rownames_to_column("area") %>%
        pivot_longer(-c(area)) %>%
        filter(area != name) %>%
        arrange(value) %>%
        mutate(
            area1 = case_when(
                area < name ~ area,
                area > name ~ name
            ),
            area2 = case_when(
                area < name ~ name,
                area > name ~ area
            )
        ) %>%
        select(-name, -area) %>%
        unique()
    
    return(sampling.points.dist)
}


# Create UMAP layout from a count matrix using Hellinger transformation
create_umap <- function(arg_count_mat){
    umap.obj <- arg_count_mat %>%
        as.matrix %>%
        t %>%
        decostand(method = "hellinger") %>%
        data.frame(check.names = F) %>%
        umap()

    umap.layout <- umap.obj$layout %>%
        data.frame()

    umap.layout$city <- rownames(umap.layout)

    return(umap.layout)
}

# Plot UMAP layout with city labels and region coloring
plot_umap  <- function(umap.df, metadata, geoCol='city', title="", max.overlaps=15, fontsize=12, palette){

    plot.data <- umap.df %>% left_join(metadata %>% select({{geoCol}}, Region) %>% distinct, by=geoCol)

    p <- ggplot(plot.data, aes(X1, X2, fill = Region)) +
        geom_point(size=3, shape=21, alpha = .9) +
        geom_text_repel(aes(label=city), max.overlaps = max.overlaps, size=fontsize) +
        scale_fill_manual(values=palette) +
        labs(x = "UMAP1", y = "UMAP2", title = title) +
        scale_x_continuous(sec.axis = sec_axis( ~.x)) +
        scale_y_continuous(sec.axis = sec_axis(~.x)) +
        guides(color=guide_legend(byrow=T, nrow=2)) + 
        theme(
            legend.position="bottom", 
            legend.text=element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5)
        )

    return(p)
}


# Run dimensionality reduction and distance analysis on ARG count data
run_dd <- function(flankodat, metadata, geoCol='city', geneCol='gene', min_copies_per_arg=1, min_args_per_site=100, verbose=F){
    if(verbose){
        print("Input arguments:")
        print(paste("geoCol:", geoCol))
        print(paste("geneCol:", geneCol))
        print(paste("Min. copies per arg:", min_copies_per_arg))
        print(paste("Min. args per site:", min_args_per_site))
    }
    
    # Create and filter count matrix
    mats <- make_mats(
        dat = flankodat, 
        geoCol = geoCol,
        geneCol = geneCol,
        min_copies_per_arg=min_copies_per_arg, 
        min_args_per_site=min_args_per_site,
        verbose = verbose
    )

    arg_count_mat <- mats$filtered

    # Run UMAP if matrix is sufficiently large
    if(all(dim(arg_count_mat)) > 5){
        umap.df <- create_umap(arg_count_mat)
        p.umap <- plot_umap(
                umap.df, 
                metadata = metadata %>% select(geoCol, Region) %>% distinct(.)
            )
        
    } else {
        umap.df <- p.umap <- NULL
    }
    
    # Compute resistance and physical distance matrices if data is non-empty
    if(all(dim(arg_count_mat)>0)){
        resDD <- create_res_dd(arg_count_mat)
        physDD <- create_physical_dd(metadata, geoCol=geoCol)
    
        # Merge resistance and physical distances
        dists <- resDD %>%
            left_join(
                physDD, 
                by = c("area1", "area2"),
                suffix=c(".res", ".m")
            ) %>%
            mutate(
                dist_m = as.numeric(value.m),
                dist_km = as.numeric(value.m) / 1000,
            )
    } else {
        resDD <- physDD <- dists <- NULL
    }

    # collect results
    res = mats
    res$umap.df <- umap.df
    res$p.umap <- p.umap
    res$dists <- dists

    return(res)

}



# Run full analysis pipeline for ResFinder and Functional ARG groups
run_all <- function(flankodat, metadata, groupCol='group', geoCol='city', geneCol='gene', min_copies_per_arg=1, min_args_per_site=100, xlim=c(NA,NA), verbose=0, title=NA){

    all_dists <- list()
    all_results <- list()

    groups <- unique(flankodat[[groupCol]])
    for (grp in groups){

        # Subset data for the current group
        grp.data <- flankodat %>% filter(group == grp)
        
        if(verbose>0){
            print(paste("Processing", grp, ":", nrow(grp.data), "x", ncol(grp.data)))
        }

        #  Run distance analyses
        grp.res <- run_dd(
            flankodat = grp.data,
            metadata = metadata,
            geoCol = geoCol,
            geneCol = geneCol,
            min_copies_per_arg = min_copies_per_arg,
            min_args_per_site = min_args_per_site,
            verbose = verbose > 1
        )
        all_results[[paste0(grp, ".res")]] <- grp.res

        # Extract and annotate distance matrix
        if(!is.null(grp.res$dists)){
            grp.dists <- grp.res$dists %>% mutate(group = grp)

            # Annotate with region data
            grp.dists <- grp.dists %>% 
                left_join(metadata %>% select(all_of(geoCol), Region) %>% distinct(), by = c("area1" = geoCol)) %>%
                left_join(metadata %>% select(all_of(geoCol), Region) %>% distinct(), by = c("area2" = geoCol), suffix=c(".1", ".2")) %>%
                na.omit()

            all_dists[[grp]] <- grp.dists

        }
    }

    # Combine all distance matrices
    combined_dists <- bind_rows(all_dists)

    # Return results
    res <- list(dists = combined_dists, results = all_results)
    return(res)
} 


# Function to extract model summary statistics
extract_model_summary <- function(model, x, z="Independent",decimals=2) {
    
    # Get adjusted R-squared
    R2 = signif(summary(model)$adj.r.squared, decimals)
    
    # Extract coefficients, filter out intercept, and format
    summary.coeffs <- as.data.frame(coefficients(summary(model)))  %>%
        rownames_to_column('term') %>%
        filter(term != '(Intercept)') %>%
        rename(Slope = Estimate, "P.value" = `Pr(>|t|)`) %>%
        mutate(
            model = z,
            R2 = R2,
            Slope = signif(Slope, decimals),
            `P.value` = signif(`P.value`, decimals),
            signif = case_when(
                P.value <= 0.0001 ~ '***',
                P.value <= 0.001 ~ '**',
                P.value <= 0.01 ~ '*',
                TRUE ~ ''
            )
        ) %>%
        mutate(Slope.str = paste0(Slope, signif)) %>%
        select(model, R2, Slope, `P.value`, signif, Slope.str)

    return(summary.coeffs)
}


# Run a Mantel test between two distance variables in a long-format data frame
run_mantel <- function(data, x, y, g1='area1', g2='area2', decimals=2, return_string=FALSE){

    # Determine matrix size based on number of unique locations
    n <- max(n_distinct(data[g1]), n_distinct(data[g2])) + 1

    # Initialize empty distance matrices
    x.mat <- y.mat <- matrix(0, n, n)

    # # Fill upper triangle of matrices with sorted distance values
    x.mat[upper.tri(x.mat)] <- data %>% arrange(!!sym(g1), !!sym(g2)) %>% pull(!!sym(x))
    y.mat[upper.tri(y.mat)] <- data %>% arrange(!!sym(g1), !!sym(g2)) %>% pull(!!sym(y))

    # Symmetrize the matrices
    x.mat <- x.mat + t(x.mat)
    y.mat <- y.mat + t(y.mat)

    # Run Mantel test
    mantel.result <- vegan::mantel(x.mat, y.mat)

    # Format the results
    results <- as.data.frame(
        t(c(`r` = signif(mantel.result$statistic, decimals), `p.value`= signif(mantel.result$signif, decimals)))
    ) %>% 
        mutate(
        signif = case_when(
                p.value <= 0.0001 ~ '***',
                p.value <= 0.001 ~ '**',
                p.value <= 0.01 ~ '*',
                TRUE ~ ''
            )
    )

    # Create a compact result string
    results.string <-paste0(results$r, results$signif)

    return(list("table" = results, "s" = results.string))
}

pretty_eq <- function(slope, intercept, r2, zz){
    
    # Determine the sign to use in the equation based on the intercept
    if(intercept < 0){sign="-"}else{sign='+'}

    # Construct the equation part with italic formatting for y
    eq = paste0("italic(y) == ", slope, "*x",sign, abs(intercept))

    
    # Construct the R-squared part with italic formatting
    r2 = paste0("italic(R)^2 == ", r2)

    
    # Combine label (zz), equation, and R-squared into a full expression string
    # This will be parsed later in ggplot2 to display formatted math
    full_eq = paste0("paste(", "'", zz, ": ', ", eq, ",', ',", r2, ")")
    return(full_eq)
}

build_lms <- function(data, x, y, z, decimals=2, palette=NA, yadj=-.1, ymin=NA, ymax=NA, xlabel=NA, ylabel=NA, title=NA, col_title=NA, calc_mantel=FALSE){
    ym <-  str_replace(y, "similarity", "dissimilarity")
    
    # Filter data
    data <- data %>%
        mutate(across(all_of(c(x, y)), ~ replace(., is.nan(.) | is.infinite(.), NA))) %>%
        filter(across(all_of(c(x, y)), ~ !is.na(.)))
    
    # Initialize an empty vector to store equations
    equations <- c()

    # Get unique values of the grouping variable z
    zs = as.vector(unique(data[[z]]))
 
    # If palette is not a list, create a default palette   
    if(!is.list(palette)){
        palette <- 2:(length(zs) + 1)
        names(palette) <- zs        
    }

    # Create the formula for the linear model    
    formula = paste(y, "~", x)

    # Fit the linear model to the entire dataset
    fittedLM <- lm(formula, data)

    # Create the equation string for the overall model
    equations[['no_scale']] <- pretty_eq(
        slope = signif(coef(fittedLM)[x], decimals),
        intercept = signif(coef(fittedLM)['(Intercept)'], decimals),
        r2 = signif(summary(fittedLM)$adj.r.squared, decimals),
        zz = 'No scale'
    )

    # Extract summary statistics for the overall model
    fittedLMs <- extract_model_summary(fittedLM, decimals=decimals)

    if(calc_mantel){
        mantelLM <- run_mantel(data = data, x = x, y = ym, decimals=decimals)
        mantel.tabular <- mantelLM$table %>% mutate(model = str_to_title('independent'))
        fittedLMs$mantel <- mantelLM$s
    }

    # Add predictions to the dataset
    data <- data %>% 
        add_predictions(fittedLM)
    
    # Create the initial ggplot object with points and a dashed line for predictions
    p <- ggplot(data, aes_string(x=x, color=z)) +
        geom_point(aes_string(y=y), alpha = .2) +
        geom_line(aes(y=pred), color="black", linewidth=1, linetype = "longdash") +
        scale_color_manual(values=palette)

    # Create a second ggplot object with smaller points
    p2 <- ggplot(data, aes_string(x=x, color=z)) +
        geom_point(aes_string(y=y), alpha = .2, size = .5) +
        geom_line(aes(y=pred), color="black", linewidth=1, linetype = "longdash") +
        scale_color_manual(values=palette)

    # Loop through each unique value of z
    for(zz in zs){
        # Filter the data for the current value of z
        subset.data <- data %>% filter(!!sym(z) == zz)

        # Fit the linear model to the subset of data
        zLM <- lm(formula, subset.data)
        
        # Get the color for the current value of z
        zz.color <- palette[zz]

        # Create the equation string for the subset model
        coeffs <- coef(zLM)

        equations[[zz]] <- pretty_eq(
            slope =  signif(coeffs[x], decimals),
            intercept =  signif(coeffs['(Intercept)'],decimals),
            r2 = signif(summary(zLM)$adj.r.squared, decimals),
            zz = zz
        )

        # Calculate 95% confidence intervals for predictions
        ci95 <- subset.data %>%
            bind_cols(
                as.data.frame(predict(zLM, subset.data, interval = "confidence", level = 0.95)) %>%
                mutate(ID=95) %>%
                remove_rownames() 
            )

        # Determine positions for annotations
        x.position <- max(min(subset.data[[x]]), min(data[[x]]))
        y.position <- mean(subset.data[[y]]) + yadj
        
        # Add predictions to the subset data
        subset.data <- subset.data %>% add_predictions(zLM, var='predZ')

        # Add lines and annotations to the first plot
        p <- p + 
            geom_line(data = subset.data, aes(y=predZ), linewidth = 1, color=darker.col(zz.color)) 
            #annotate("text", x = x.position, y = y.position, label = eq, hjust = 0, vjust = 1, color = darker.col(zz.color, 75))
        
        # Add lines, annotations, and confidence intervals to the second plot
        p2 <- p2 + 
            geom_line(data = subset.data, aes(y=predZ), linewidth = 1, color=darker.col(zz.color)) + 
            #annotate("text", x = x.position, y = y.position + yadj*3, label = eq2, hjust = 0, vjust = 1, color = darker.col(zz.color, 75)) +
            geom_ribbon(data=ci95, aes(ymin=lwr, ymax=upr), fill = darker.col(zz.color, 20), alpha = .4, linewidth = 0)

        
        fittedLMz <- extract_model_summary(zLM, decimals=decimals, z = str_to_title(zz))

        if(calc_mantel){
            mantelLM <- run_mantel(data = subset.data, x = x, y = str_replace(y, "similarity", "dissimilarity"), decimals=decimals)
            mantel.tabular <- rbind(mantel.tabular, mantelLM$table %>% mutate(model = str_to_title(zz)))
            fittedLMz$mantel <- mantelLM$s
        }
        fittedLMs <- rbind(
            fittedLMs, fittedLMz
        )
    }

    # Fit a linear model with interaction terms
    if(length(zs)>1){
        fittedLMz <- lm(paste0(y, "~", x, "*", z), data)

        # Get pairwise comparisons of trends
        lmz.lst <- emmeans::emtrends(fittedLMz, z, var=x)
        pairs.df <- as.data.frame(pairs(lmz.lst, adjust = "tukey"))

        # Create a bar plot for pairwise comparisons
        p.bar <- pairs.df %>%
            #mutate(contrast = str_wrap(contrast, width = 10)) %>%
            mutate(contrast = str_replace_all(contrast, " ", "\n")) %>%
            ggplot(aes( x = contrast, y = estimate, fill = p.value < 0.05)) +
            geom_bar(stat = "identity") +
            geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.2) +
            labs(title = "Pairwise Comparisons of Trends", x = "Comparison", y = "Estimated Difference") +
            theme_minimal() +
            scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = paste("P", "<", 0.05)) + 
            theme(
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5)
        )
    } else {
        p.bar <- ggplot() + theme_void()
        pairs.df <- NA
    }
    
    # Tidy up the plot titles and labels
    if(is.na(col_title)){
        col_title <- str_replace(string = z, pattern = "_", replacement = " ") %>% str_to_title()
    }
    
    p <- p +
        guides(
            col=guide_legend(
                col_title, 
                override.aes = list(shape=16, alpha=1)
            )
        )  + 
        theme(
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)
        ) +
        labs(title=title)

    p2 <- p2 +
        guides(
            col=guide_legend(
                col_title, 
                override.aes = list(shape=16, alpha=1)
            )
        )  + 
        theme(
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)
        ) +
        labs(title=title)

    # Add axis labels and limits if specified
    if(!is.na(ymin)){
        p <- p + ylim(ymin, ymax)
        p2 <- p2 + ylim(ymin, ymax)
    }

    if(!is.na(xlabel)){
        p <- p + labs(x=xlabel)
        p2 <- p2 + labs(x=xlabel)

    }

    if(!is.na(ylabel)){
        p <- p + labs(y=ylabel)
        p2 <- p2 + labs(y=ylabel)

    }  

    if(calc_mantel){
        fittedLMs.tab <- fittedLMs %>% 
        column_to_rownames('model') %>% 
        select(R2, Slope.str, mantel) %>%
        rename(Slope = Slope.str, Mantel = mantel)
    } else {
        fittedLMs.tab <- fittedLMs %>% 
        column_to_rownames('model') %>% 
        select(R2, Slope.str) %>%
        rename(Slope = Slope.str)
    }
    fittedLMs.tabular <- gridExtra::tableGrob(fittedLMs.tab, theme = gridExtra::ttheme_minimal())

    # Combine the second plot with the bar plot
    p.bottom <- ggarrange(fittedLMs.tabular, p.bar, ncol = 2, nrow = 1)
    #p.combined <- ggarrange(p2, p.bar, ncol=1, heights=c(2, 1))
    p.combined <- ggarrange(p2, p.bottom, ncol=1, heights=c(2, 1))

    all.results <- 
        list(
            "p" = p, 
            "p2" = p.combined, 
            "scale.independent.model" = fittedLM, 
            "scale.dependent.model" = fittedLMz, 
            "pairs" = pairs.df, 
            "fit.summaries" = fittedLMs,
            "equations" = equations
        )
    if(calc_mantel){
        all.results$mantel.results <- mantel.tabular
        }
    return(all.results)
}

add_plot_equations <- function(p, eqns, xmin, ymin, jitter=10, size = 3){

    n = length(eqns)
    for(i in 1:n){
        y_pos = ymin + (i-1)/jitter
        p <- p + annotate(
            "text",
            x = xmin,
            y = y_pos,
            label = rev(eqns)[[i]],
            parse = T,
            hjust = 0,
            size = size
        )
    }
    return(p)
}
