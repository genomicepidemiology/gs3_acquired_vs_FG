library(dplyr)
library(tidyr)


# Calculate the distribution of ARGs across geographic areas
process_spread <- function(data, geo_level, arg_level, group, res_class = 'All'){

    # Filter by group and resistance class
    sel.data <- data[which(data$group == group),]

    if(res_class != 'All'){
        sel.data <- sel.data[which(sel.data$class == res_class),]
    }
    if(any(dim(sel.data) == 0)){
        return(c())
    }
    # Count presence of ARGs across geographic units
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

# Generate spread statistics across multiple groupings and resistance classes
generate_spread <- function(dat1, dat2, geo_levels, resClasses, groups=c("ResFinder", "Functional"), arg_levels=c("gene", "variant")){

    spreads <- c()
    for(group in groups){
        for(arg_level in arg_levels){
            for(geo_level in geo_levels){
                print(paste(group, arg_level, geo_level))

                # Process 'All' resistance classes from dat1
                spread.data <- process_spread(
                    data = dat1,
                    geo_level = geo_level,
                    arg_level = arg_level,
                    group = group
                )

                spreads <- rbind(spreads, spread.data)

                # Process specific resistance classes from dat2
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

    # Format factor levels for plotting or analysis
    spreads <- spreads %>%
        mutate(
            group = factor(group, level=groups),
            geo_level = factor(geo_level, levels=geo_levels),
            arg_level = factor(arg_level, levels = arg_levels),
            res_class = factor(res_class, levels = c("All", resClasses))
        )

    return(spreads)
}