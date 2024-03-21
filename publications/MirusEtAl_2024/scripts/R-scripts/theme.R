library(ggplot2)

custom_colours <- list(
	lit_colours = c("#12becc", "#f58080", "#0059ec", "#fcab79", "#6a00e5"),
	lit_colours_all = c("#12becc", "#f58080", "#0059ec", "#ffd266", "#6a00e5", "#fcab79") 
)

custom_palettes = function(name, n, all_palettes = custom_colours, type = c("discrete", "continuous")) {
	palette = all_palettes[[name]]
  	if (missing(n)) {
		n = length(palette)
    	}
    	type = match.arg(type)
    	out = switch(type, continuous = grDevices::colorRampPalette(palette)(n), discrete = palette[1:n])
      	structure(out, name = name, class = "palette")
}

scale_colour_custom_d <- function(name) {
	ggplot2::scale_colour_manual(values = custom_palettes(name, type = "discrete"))
}

scale_colour_custom_c <- function(name) {
	ggplot2::scale_colour_gradientn(colours = custom_palettes(name, type = "continuous"))
}

scale_fill_custom_d <- function(name) {
	ggplot2::scale_fill_manual(values = custom_palettes(name, type = "discrete"))
}

scale_fill_custom_c <- function(name) {
	ggplot2::scale_fill_gradientn(colours = custom_palettes(name, type = "continuous"))
}

custom_theme <- theme_bw() +
	theme(
  		#plot.title = element_blank(),
  		plot.subtitle = element_blank(),
  		axis.title = element_text(size = 24),
  		axis.text = element_text(size = 20),
  		legend.title = element_text(size = 24),
  		legend.text = element_text(size = 20),
		strip.text.x = element_text(size = 24),
		strip.text.y = element_text(size = 24, margin = margin(0.2, 0.2, 0.2, 0.2, "cm")),
		panel.grid = element_blank(),
		strip.background = element_rect(fill = "white", colour = "grey")
	)

custom_theme_large <- theme_bw() +
	theme(
  		#plot.title = element_blank(),
  		plot.subtitle = element_blank(),
  		axis.title = element_text(size = 36),
  		axis.text = element_text(size = 32),
  		legend.title = element_text(size = 36),
  		legend.text = element_text(size = 32),
		strip.text.x = element_text(size = 32),
		strip.text.y = element_text(size = 32, margin = margin(0.2, 0.2, 0.2, 0.2, "cm")),
		panel.grid = element_blank(),
		strip.background = element_rect(fill = "white", colour = "grey")
	)
