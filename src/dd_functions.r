make_count_matrix <- function(dat, geoCol='city', geneCol='gene'){
    mat <- dat %>%
        select(all_of(geneCol), all_of(geoCol)) %>%
        table %>%
        as.data.frame.matrix
    return(mat)
}

filter_count_matrix <- function(mat, min_copies_per_arg=1, min_args_per_site=100) {
    return(mat[rowSums(mat)>min_copies_per_arg, colSums(mat)>min_args_per_site])
}

make_mats <- function(dat, geoCol='city', geneCol='gene', min_copies_per_arg=1, min_args_per_site=100, verbose=F){
    unfiltered.mat <- make_count_matrix(
        dat = dat, 
        geoCol = geoCol,
        geneCol= geneCol
    )
    
    filtered.mat <- filter_count_matrix(
        mat= unfiltered.mat, 
        min_copies_per_arg=min_copies_per_arg, 
        min_args_per_site=min_args_per_site
    )

    if(verbose){
        print(paste("Dimension of input data:", nrow(dat), ncol(dat)))
        print(paste("Dimension of unfiltered count matrix:", nrow(unfiltered.mat), ncol(unfiltered.mat)))
        print(paste("Dimension of filtered count matrix:", nrow(filtered.mat), ncol(filtered.mat)))
        
    }

    return(list("mat" = unfiltered.mat, "filtered" = filtered.mat))
}

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

create_physical_dd <- function(metadata, geoCol){

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

    sampling.points.dist <- st_distance(sampling.points, sampling.points) %>%
        data.frame

    colnames(sampling.points.dist) <- rownames(sampling.points.dist) <- sampling.points[[geoCol]]

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

run_dd <- function(flankodat, metadata, geoCol='city', geneCol='gene', min_copies_per_arg=1, min_args_per_site=100, verbose=F){
    if(verbose){
        print("Input arguments:")
        print(paste("geoCol:", geoCol))
        print(paste("geneCol:", geneCol))
        print(paste("Min. copies per arg:", min_copies_per_arg))
        print(paste("Min. args per site:", min_args_per_site))
    }
    
    mats <- make_mats(
        dat = flankodat, 
        geoCol = geoCol,
        geneCol = geneCol,
        min_copies_per_arg=min_copies_per_arg, 
        min_args_per_site=min_args_per_site,
        verbose = verbose
    )

    arg_count_mat <- mats$filtered
    if(all(dim(arg_count_mat)) > 5){
        umap.df <- create_umap(arg_count_mat)
        p.umap <- plot_umap(
                umap.df, 
                metadata = metadata %>% select(geoCol, Region) %>% distinct(.)
            )
        
    } else {
        umap.df <- p.umap <- NULL
    }
    
    if(all(dim(arg_count_mat)>0)){
        resDD <- create_res_dd(arg_count_mat)
        physDD <- create_physical_dd(metadata, geoCol=geoCol)
    
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

    res = mats
    res$umap.df <- umap.df
    res$p.umap <- p.umap
    res$dists <- dists

    return(res)

}

run_all <- function(flankodat, metadata, geoCol='city', geneCol='gene', min_copies_per_arg=1, min_args_per_site=100, xlim=c(NA,NA), verbose=0, title=NA){

    rf.data <- flankodat %>% filter(group == 'ResFinder')
    if(verbose > 0){
        print(paste("ResFinder data input:", nrow(rf.data), ncol(rf.data)))
    }
    rf.res <- run_dd(
        flankodat = rf.data,
        metadata = metadata,
        geoCol = geoCol,
        geneCol = geneCol,
        min_copies_per_arg = min_copies_per_arg,
        min_args_per_site = min_args_per_site,
        verbose=verbose > 1
    )
    func.data <- flankodat %>% filter(group == 'Functional')
    if(verbose > 0){print(paste("Functional data input:", nrow(func.data), ncol(func.data)))}
    func.res <- run_dd(
        flankodat = flankodat %>% filter(group == 'Functional'),
        metadata = metadata,
        geoCol = geoCol,
        geneCol = geneCol,
        min_copies_per_arg = min_copies_per_arg,
        min_args_per_site = min_args_per_site
    )

    res.dists <- rf.res$dists 
    func.dists <- func.res$dists
    if(verbose > 0){
        print(paste("ResFinder distance matrix:", nrow(res.dists), ncol(res.dists)))
        print(paste("Functional distance matrix:", nrow(func.dists), ncol(func.dists)))
    }

    if(!is.null(res.dists)){res.dists <- res.dists %>% mutate(group = 'ResFinder') }
    if(!is.null(func.dists)){func.dists <- func.dists %>% mutate(group = 'Functional') }
    
    dists <- rbind(res.dists, func.dists) 
    if(!is.null(dists)){
        dists <- dists %>%
        left_join(
            metadata %>% select(all_of(geoCol), Region) %>% distinct(.),
            by = c("area1" = geoCol)
        ) %>%
        left_join(
            metadata %>% select(all_of(geoCol), Region) %>% distinct(.),
            by = c("area2" = geoCol),
            suffix=c(".1", ".2")
        ) %>%
        na.omit

    }
    
    allRes <- c()

    allRes$dists  <- dists
    allRes$rf.res <- rf.res
    allRes$func.res <- func.res
    return(allRes)   
}

darker.col = function(color, how.much = 30){
  colorRampPalette(c(color, "black"))(100)[how.much]
}

simulate_data <- function(min_distance, max_distance, model){
    # Create the simulated datasets
    sim_data_FALSE <- data.frame(distance = seq(min_distance, max_distance, length.out = 100))
    sim_data_TRUE <- data.frame(distance = seq(min_distance, max_distance, length.out = 100))
    
    # Add the within column
    sim_data_FALSE$within <- FALSE
    sim_data_TRUE$within <- TRUE
    
    # Calculate the log of distances
    sim_data_FALSE$log_distance <- log(sim_data_FALSE$distance)
    sim_data_TRUE$log_distance <- log(sim_data_TRUE$distance)
    
    # Use the model coefficients to predict log_dissimilarity
    coefficients <- coefficients(model)
    # For within = FALSE
    intercept <- coefficients[["(Intercept)"]] #-0.47875
    slope_log_distance <- coefficients[["log_distance"]]# 0.00650
    
    sim_data_FALSE$log_dissimilarity <- intercept + slope_log_distance * sim_data_FALSE$log_distance
    
    # For within = TRUE
    intercept_withinTRUE <- coefficients[["withinTRUE"]] #-0.65397
    slope_interaction <- coefficients[["log_distance:withinTRUE"]] #0.06725
    
    sim_data_TRUE$log_dissimilarity <- intercept + intercept_withinTRUE + 
                                       (slope_log_distance + slope_interaction) * sim_data_TRUE$log_distance
    
    # Combine the datasets if needed
    sim_data <- rbind(sim_data_FALSE, sim_data_TRUE)
    
    # View the simulated data
    return(sim_data)

}

fit_dd_model <- function(data, x='dist_km', y='value.res', ymin=NA, ymax=NA, title="", col_title="Within-Region", add_sim=F, xlabel=NA, ylabel=NA, palette=NA){

    if(!is.list(palette)){
        palette <- list("TRUE" = "blue", "FALSE" = "red")
    }

    data <- data %>% 
        filter(!is.na(!!sym(x)), !is.na(!!sym(y))) %>%
        mutate(
            distance = !!sym(x),
            dissimilarity = !!sym(y),
            within = Region.1 == Region.2,
            log_distance = log(distance),
            log_dissimilarity = log(dissimilarity)
        ) %>%
        mutate(across(everything(), ~ replace(., is.nan(.) | is.infinite(.), NA))) %>%
        filter(across(everything(), ~ !is.na(.)))

    n <- max(c(n_distinct(data$area1), n_distinct(data$area2))) + 1
    dist.mat <- matrix(0, n, n)
    diss.mat <- matrix(0, n, n)
    
    dist.mat[upper.tri(dist.mat)] <- data %>% arrange(area1, area2) %>% pull(dist_km)
    diss.mat[upper.tri(diss.mat)] <- data %>% arrange(area1, area2) %>% pull(dissimilarity)

    dist.mat <- dist.mat + t(dist.mat)
    diss.mat <- diss.mat + t(diss.mat)

    mantel_results <- mantel(apply(dist.mat, 2, rank), apply(diss.mat, 2, rank))
    mantel_string <- paste("Mantel r =", round(mantel_results$statistic,4), "(p-value =", format(mantel_results$signif, scientific=T), ")")

    # model
    model.all <- lm(log_dissimilarity ~ log_distance, data=data)
    data$predicted <- predict(model.all, type = "response")

    # regionality models
    model.region <- lm(log_dissimilarity ~ log_distance*within, data = data)
    data$predicted_region <- predict(model.region, type = "response")

    # get the coefficients for the two slopes
    coefficients <- summary(model.region)$coefficients

    slope_within_FALSE = coefficients["log_distance", "Estimate"]
    pval_within_FALSE  = coefficients["log_distance", "Pr(>|t|)"]
    slope_within_TRUE = coefficients["log_distance", "Estimate"] + coefficients["log_distance:withinTRUE", "Estimate"]
    pval_within_TRUE = coefficients["log_distance:withinTRUE", "Pr(>|t|)"]

    string_within_FALSE = paste(
        "Slope for between:", format(round(slope_within_FALSE, 4), scientific=T), 
        "(p-value =", format(pval_within_FALSE, scientific=T), ")"
    )
    string_within_TRUE = paste(
        "Slope for within:", format(round(slope_within_TRUE, 4), scientific = T), 
        "(p-value =", format(pval_within_TRUE, scientific=T), ")"
    )

    #return(sim.data)
    p <- ggplot(data, aes(x=log_distance)) +
        geom_point(aes(y=log_dissimilarity, color=within), alpha = .1, shape=1) +
        scale_color_manual(values=palette) +
        geom_line(aes(y=predicted), size=1, linetype = "dashed") +
        geom_line(
            data=data %>% filter(within == T), 
            aes(y=predicted_region),
            color=darker.col(palette$`TRUE`), 
            size=1)+
        geom_line(
            data=data %>% filter(within == F), 
            aes(y=predicted_region),
            color=darker.col(palette$`FALSE`), 
            size=1
        ) +
        labs(
            title=title
        ) + 
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
        annotate("text", x = -Inf+0.2, y=-Inf, hjust = 0.001, vjust = -2.2, label = mantel_string) +
        annotate("text", x = -Inf+0.2, y=-Inf, hjust = 0.001, vjust = -1.2, label = string_within_TRUE, color=darker.col(palette$`TRUE`, 50))+ 
        annotate("text", x = -Inf+0.2, y=-Inf, hjust = 0.001, vjust = -.1, label = string_within_FALSE, color=darker.col(palette$`FALSE`, 50))
    
    if(!is.na(ymin)){
        p <- p + ylim(ymin, ymax)
    }

    if(!is.na(xlabel)){
        p <- p + labs(x=xlabel)
    } else {
        p <- p + labs(x="Ln(Geodesic distance) (km)")
    }
    if(!is.na(ylabel)){
        p <- p + labs(y=ylabel)
    } else {
        p <- p + labs(y="Ln(Resistome Dissimilarity)")
    }

    if(add_sim){
        sim.data <- simulate_data(
            min_distance = min(data$distance), 
            max_distance=max(data$distance), 
            model=model.region
        )

        p <- p +
            geom_line(
                data = sim.data %>% filter(within == F),
                mapping = aes(y=log_dissimilarity),
                color=darker.col(palette$`FALSE`),
                size = 1,
                alpha = .5,
                linetype = "dashed"
            ) +
            geom_line(
                data = sim.data %>% filter(within == T),
                mapping = aes(y=log_dissimilarity),
                color=darker.col(palette$`TRUE`),
                size = 1,
                alpha = .5,
                linetype = "dashed"
            )
    }

    return(list("p" = p, "model.all" = model.all, "mantel" = mantel_results, "model.region" = model.region))
}

make_clr_dist <- function(mat, id_col='genepid',method='Euclidean', drop_columns=c()){

    dist.mat <- as.matrix(dist(mat %>% select(-drop_columns) %>% column_to_rownames(id_col), method))

    dist.df <- as_tibble(dist.mat, rownames=NA) %>%
        rownames_to_column(id_col) %>%
        pivot_longer(-c(id_col), values_to = paste0('dist.', method), names_to = paste0(id_col, '_2'))
    return(dist.df)
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

fit_dd_model2 <- function(data, x='dist_km', y='value.res', dissim=TRUE, ymin=NA, ymax=NA, title="", col_title="Spatial scale", add_sim=F, xlabel=NA, ylabel=NA, palette=NA){

    if(!is.list(palette)){
        palette <- list("Within Country" = "#F8766D", "Within Region" = "#00BFC4", "Between Regions" = "#C77CFF")
    }

    if(dissim){y_col = 'log_dissimilarity'}else{y_col='log_similarity'}
    
    data <- data %>%
        mutate(
            distance = !!sym(x),
            dissimilarity = !!sym(y),
            similarity = 1-!!sym(y),
            log_distance = log(distance),
            log_dissimilarity = log(dissimilarity),
            log_similarity = log(similarity),
            spatial_scale = factor(
                    case_when(
                    country.1 == country.2 ~ 'Within Country',
                    country.1 != country.2 & Region.1 == Region.2 ~ 'Within Region',
                    TRUE ~ 'Between Regions'
                ),
                levels=c('Within Country', 'Within Region', 'Between Regions')
            )
        ) %>%
        mutate(across(everything(), ~ replace(., is.nan(.) | is.infinite(.), NA))) %>%
        filter(across(everything(), ~ !is.na(.))) %>%
        arrange(desc(spatial_scale))

    # subset spatial_scales to those in legend?
    data <- data[which(data$spatial_scale %in% names(palette)),]

    n <- max(c(n_distinct(data$area1), n_distinct(data$area2))) + 1
    dist.mat <- matrix(0, n, n)
    diss.mat <- matrix(0, n, n)

    dist.mat[upper.tri(dist.mat)] <- data %>% arrange(area1, area2) %>% pull(dist_km)
    if(dissim){
        diss.mat[upper.tri(diss.mat)] <- data %>% arrange(area1, area2) %>% pull(dissimilarity)
    } else {
        diss.mat[upper.tri(diss.mat)] <- data %>% arrange(area1, area2) %>% pull(dissimilarity)

    }

    dist.mat <- dist.mat + t(dist.mat)
    diss.mat <- diss.mat + t(diss.mat)

    mantel_results <- mantel(apply(dist.mat, 2, rank), apply(diss.mat, 2, rank))
    mantel_string <- paste("Mantel r =", round(mantel_results$statistic,4), "(p-value =", format(mantel_results$signif, scientific=T), ")")

    # model
    if(dissim){
        formula = "log_dissimilarity ~ log_distance"
    } else {
        formula = "log_similarity ~ log_distance"
    }
    model.all <- lm(formula, data=data, na.action=na.exclude)
    data$predicted <- predict(model.all, type = "response")

    pred_data <- c()
    spatial.models <- list("all" = model.all)
    for(spatial_group in names(palette)){
        spatial_data <- data %>% filter(spatial_scale == spatial_group)
        if(nrow(spatial_data) > 0){
            spatial_model <- lm(formula, data=spatial_data,na.action=na.exclude)
            spatial_data$predicted <- predict(spatial_model, type="response")
    
            pred_data <- rbind(pred_data, spatial_data)
            spatial.models[[spatial_group]] <- spatial_model
        }
    }
    #data$predicted_region <- predict(model.region, type = "response")


    p <- ggplot(data, aes(x=log_distance)) +
        geom_point(aes(y=!!sym(y_col), color=spatial_scale), alpha = .3, shape = 1) +
        scale_color_manual(values=palette)+
        geom_line(aes(y=predicted), size=1, linetype = "dashed") +
        geom_line(
            data=pred_data %>% filter(spatial_scale == 'Within Country'), 
            aes(y=predicted),
            color = darker.col(palette$`Within Country`),
            size = 1
        )+
        geom_line(
            data=pred_data %>% filter(spatial_scale == 'Between Regions'), 
            aes(y=predicted),
            color = darker.col(palette$`Between Regions`),
            size = 1
        )+
        geom_line(
            data=pred_data %>% filter(spatial_scale == 'Within Region'), 
            aes(y=predicted),
            color = darker.col(palette$`Within Region`),
            size = 1
        ) + 
        guides(
            col=guide_legend(
                col_title, 
                override.aes = list(shape=16, alpha=1)
            )
        ) + 
        theme(
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)
        )+
        labs(
            title=title
        )

    if(!is.na(ymin)){
        p <- p + ylim(ymin, ymax)
    }

    if(!is.na(xlabel)){
        p <- p + labs(x=xlabel)
    } else {
        p <- p + labs(x="Ln(Geodesic distance) (km)")
    }
    if(!is.na(ylabel)){
        p <- p + labs(y=ylabel)
    } else {
        if(dissim){
            p <- p + labs(y="Ln(Resistome Dissimilarity)")
        } else {
            p <- p + labs(y="Ln(Resistome Similarity)")            
        }
    }
    
    return(list("p" = p, "models" = spatial.models, "mantel" = mantel_results, "data" = data))

}

summarise_dd_model_slopes <- function(models, decimals=4, col_name='vv'){
    model.coefficients <- c()
    for(model_name in names(models)){
        coefficients <- as.data.frame(coefficients(summary(models[[model_name]])))%>%
            rownames_to_column('coefficient') %>%
            filter(coefficient != '(Intercept)') %>%
            mutate(
                spatial_scale = model_name,
                star = case_when(
                    `Pr(>|t|)` < 0.0001 ~ '***',
                    `Pr(>|t|)` < 0.001 ~ '**',
                    `Pr(>|t|)` < 0.01 ~ '*',
                    `Pr(>|t|)` < 0.05 ~ '.',
                    TRUE ~ ' '
                ),
                vv = paste0(
                    format(round(Estimate, decimals), scientific=F),
                    star,
                    ' (p-value = ', 
                    format(`Pr(>|t|)`, scientific=T),
                    ')'
                    )

            ) %>%
            select(spatial_scale, vv)

            

        model.coefficients <- rbind(model.coefficients, coefficients)
    }
    colnames(model.coefficients)[2] <- col_name

    return(model.coefficients)
}

run_mantel <- function(mat, metadata, method='bray', transform = F){


    if(transform){
        mat <- t(mat)
        mat <- mat %>% decostand(method="hellinger")
    } 
    res.dists <- proxy::dist(mat, method = method)

    cities <- metadata %>%
        mutate(
            city = case_when(
                city == "Santiago" & country == "Spain" ~ "Santiago_ESP",
                TRUE ~ city)
        ) %>%
        select(city, lon, lat) %>%
        unique() %>% na.omit %>%
        group_by(city) %>%
        summarise(lon = median(lon, na.rm = T), lat = median(lat, na.rm = T)) %>%
        st_as_sf(coords = c("lon", "lat"), crs = "WGS84") %>%
        filter(city %in% rownames(mat))

    phys.dists <- st_distance(cities, cities)
    colnames(phys.dists) <- rownames(phys.dists) <- cities$city    
    return(mantel(phys.dists, res.dists))
    
}