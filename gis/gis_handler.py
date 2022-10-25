import rasterio
import rasterio.plot
import matplotlib.pyplot as plt

color_map = 'turbo'

class Handler:
    
    def __init__(self, files):
        self.raster_data = {}
        load_raster_files(files)     
    
    def load_raster_files(self, files):
        for key, src in files:
            self.raster_data[key] = Raster(key, src)
    
    def list_rasters(self):
        return list(self.raster_data.keys())
    
    def display_raster(self, key):
        raster = self.raster_data[key]
        fig, ax = plt.subplots()
        ax = rasterio.plot.show(raster.data, extent=raster.extent, ax=ax, cmap=color_map)
    
    def
        
class Raster:
    
    def __init__(self, name, src):
        self.name = name
        self.src = src
        
        self.data = rasterio.open(self.src)
        self.extent = [self.data.bounds[0], self.data.bounds[2], self.data.bounds[1], self.data.bounds[3]]