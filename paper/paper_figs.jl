using QGDipoles, Plots

grid = CreateGrid(512, 512, 10, 10)

# First panel (LCD):

ψ1, = CreateLCD(grid, 1, 1)

heatmap(grid.x, grid.y, transpose(ψ1);
	colormap = :balance,
	guidefontsize=20,
	tickfontsize=16,
	aspect_ratio=1,
	xlims = (-5,5),
	ylims = (-5,5),
	xlabel = "x",
	ylabel = "y",
	size = (600,400))

savefig("Fig1a.svg")

# Second Panel (LRD):

ψ2, = CreateLRD(grid, 1, 1, 1, 1)

heatmap(grid.x, grid.y, transpose(ψ2);
	colormap = :balance,
	guidefontsize=20,
	tickfontsize=16,
	aspect_ratio=1,
	xlims = (-5,5),
	ylims = (-5,5),
	xlabel = "x",
	ylabel = "y",
	size = (600,400))

savefig("Fig1b.svg")

# Third Panel (Mode 2 SQG):

ψ3, = CreateModonSQG(grid, 12; K₀=8)

heatmap(grid.x, grid.y, transpose(ψ3);
	colormap = :balance,
	guidefontsize=20,
	tickfontsize=16,
	aspect_ratio=1,
	xlims = (-5,5),
	ylims = (-5,5),
	xlabel = "x",
	ylabel = "y",
	size = (600,400))

savefig("Fig1c.svg")

nothing
