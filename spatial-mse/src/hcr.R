hcr <- function(B, l = 19958, u = 38555, hmin = 0, hmax = .2) {
    h <- ifelse(B < l, hmin, 
        ifelse(B >= u, hmax,
            hmax/(u-l)*(B-l) 
        ) 
    )
    return(h)
}