from shiny import App, ui, render
import plotly.express as px
import pandas as pd
import numpy as np
from shinywidgets import output_widget, render_plotly

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_slider(
            "m",
            "Natural Mortality",
            min=0.15,
            max=0.5,
            value=0.25,
        ),
    ),
    ui.layout_columns(
        ui.value_box(
            "Interesting bit", ui.output_ui("interesting_bit"),
        ),
        ui.value_box(
            "Fascinating result", ui.output_ui("fascinating_result"), 
        ),
        ui.value_box(
            "Fabulous insight",
            ui.output_ui("fabulous_insight"),
        ),
        fill=False,
    ),
    ui.layout_columns(
        ui.card(
            ui.card_header("Comparison to base model biomass"),
            output_widget("compare_biomass"),
            full_screen=True,
        ),
    ),
)

def server(input, output, session):
    @render.text
    def greeting():
        return f"Hello, {input.name()}!"
    
    @render_plotly
    def compare_biomass():
        # Sample data
        df = {
            "year": range(1, 50),
            "error": np.random.uniform(-5, 5, 49)
        }        
        fig = px.line(df, x="year", y="error")
        return fig

app = App(app_ui, server)