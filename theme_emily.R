#theme_emily <- function(grid=TRUE, legend.position=NA, base_family='Lato Light', highlight_family='Lato') {
theme_emily <- function(grid=TRUE, legend.position=NA, base_family='Arial', highlight_family='Arial') {
    
  # credit to Henrik Lindberg https://github.com/halhen
  th <- ggplot2::theme_minimal(base_family = base_family, base_size = 12)
  
  th <- th + theme(text = element_text(color='#333333'))
  
  th <- th + theme(legend.background = element_blank())
  th <- th + theme(legend.key = element_blank())
  
  # Straight out of hrbrthemes
  if (inherits(grid, "character") | grid == TRUE) {
    # th <- th + theme(panel.grid=element_line(color="#cccccc", size=0.3))
    # th <- th + theme(panel.grid.major=element_line(color="#eaeaea", size=0.3))
    # th <- th + theme(panel.grid.minor=element_line(color="#eaeaea", size=0.15))
    # 
    th <- th + theme(panel.grid=element_line(color="#cccccc", size=0.3))
    th <- th + theme(panel.grid.major=element_blank())
    th <- th + theme(panel.grid.minor=element_blank())
    
    th <- th + theme(axis.line=element_line(color="grey50", size=0.4))
   # th <- th + theme(axis.line=element_line(color="#eaeaea", size=0.7))
    
    if (inherits(grid, "character")) {
      if (regexpr("X", grid)[1] < 0) th <- th + theme(panel.grid.major.x=element_blank())
      if (regexpr("Y", grid)[1] < 0) th <- th + theme(panel.grid.major.y=element_blank())
      if (regexpr("x", grid)[1] < 0) th <- th + theme(panel.grid.minor.x=element_blank())
      if (regexpr("y", grid)[1] < 0) th <- th + theme(panel.grid.minor.y=element_blank())
    }
    
  } else {
    th <- th + theme(panel.grid=element_blank())
  }
  
  th <- th + theme(axis.text = element_text(family=highlight_family))
  th <- th + theme(strip.text = element_text(family=highlight_family, face = "bold"))
  #th <- th + theme(axis.ticks = element_blank())
  
  th <- th + theme(axis.text.x=element_text(margin=margin(t=5)))
  th <- th + theme(axis.text.y=element_text(margin=margin(r=5)))
  
  # th <- th + theme(plot.title = element_text(family="Playfair Display"),
  #      plot.subtitle = element_text(margin=margin(b=15), family="Playfair Display"),
  #      plot.caption = element_text(face='italic', size=10))
  
  if (!is.na(legend.position)) th <- th + theme(legend.position = legend.position)
  
  return (th)
  
}


