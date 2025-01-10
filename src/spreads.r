library(dplyr)
library(tidyr)

process_spread <- function(data, geo_level, arg_level, group, res_class = 'All'){

    sel.data <- data[which(data$group == group),]

    if(res_class != 'All'){
        sel.data <- sel.data[which(sel.data$class == res_class),]
    }
    if(any(dim(sel.data) == 0)){
        return(c())
    }
    spread <- sel.data %>% 
        select({{geo_level}}, {{arg_level}}) %>%
        table %>%
        as_tibble() %>%
        mutate(present = (n>0)*1) %>%
        select(sequence = 2, present) %>%
        group_by(sequence) %>%
        summarise(areaCount = sum(present)) %>%
        arrange(-areaCount) %>%
        pull(areaCount) %>%
        table %>%
        as_tibble() %>%
        rename(numberAreas = 1) %>%
        filter(numberAreas > 0) %>%
        mutate(
            numberAreas = as.numeric(numberAreas),
            geo_level = geo_level,
            arg_level = arg_level,
            group = group,
            res_class = res_class
        )

    return(spread)
}

generate_spread <- function(dat1, dat2, geo_levels, resClasses, groups=c("ResFinder", "Functional"), arg_levels=c("gene", "variant")){

    spreads <- c()
    for(group in groups){
        for(arg_level in arg_levels){
            for(geo_level in geo_levels){
                print(paste(group, arg_level, geo_level))

                spread.data <- process_spread(
                    data = dat1,
                    geo_level = geo_level,
                    arg_level = arg_level,
                    group = group
                )

                spreads <- rbind(spreads, spread.data)

                for(resClass in resClasses){
                    spread.data <- process_spread(
                        data = dat2,
                        geo_level = geo_level,
                        arg_level = arg_level,
                        group = group,
                        res_class = resClass
                    )
                    
                    spreads <- rbind(spreads, spread.data)

                }
                
            }
        }
    }

    spreads <- spreads %>%
        mutate(
            group = factor(group, level=groups),
            geo_level = factor(geo_level, levels=geo_levels),
            arg_level = factor(arg_level, levels = arg_levels),
            res_class = factor(res_class, levels = c("All", resClasses))
        )

    return(spreads)
}

plot_spread <- function(spreadData, palette, x='numberAreas', y='n', fill='group', wrap="geo_level ~ arg_level", dodge_width = .5, title="", marker_size=1){

    wrap_labels <- as.vector(str_split(wrap, " ~ ", simplify = T))

    # Create a new data frame with dodged x positions
    spreadData_dodged <- spreadData %>%
        group_by(!!sym(fill)) %>%
        mutate(
            dodged_x = as.numeric(as.factor(!!sym(x))) + (as.numeric(as.factor(!!sym(fill))) - 1) * dodge_width
        ) %>%
        mutate(across(all_of(wrap_labels), str_to_title))
    
    dodge <- position_dodge(width = dodge_width)

    p <- ggplot(spreadData_dodged, aes(x=dodged_x, y=!!sym(y), color=!!sym(fill))) +
        geom_segment(aes(x=dodged_x, end=dodged_x, y=0, yend=!!sym(y)), alpha = .3) +
        geom_point(size=marker_size, alpha = .75, position = dodge) +
        scale_color_manual(values = palette) +
        facet_wrap(wrap, scales="free", ncol=2) +
        scale_y_log10() +
        theme_minimal() +
        labs(
            y="log(Number of sequences with prevalence)",
             x = "Number of Areas"
        ) +
        ggtitle(title)

    return(p)
}