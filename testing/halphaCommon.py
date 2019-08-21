#from PyQt5 import  QtWidgets
from PyQt5 import QtCore,QtWidgets, QtGui
#from PyQt5.Qtcore import  Qt
#from ginga.qtw.QtHelp import QtGui, QtCore

from halphav3 import Ui_MainWindow
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
from ginga import colors
from ginga.canvas.CanvasObject import get_canvas_types
from ginga.misc import log
from ginga.util.loader import load_data
from astropy.wcs import WCS

class cutout_image():
    #key_pressed = QtCore.pyqtSignal(str)
    def __init__(self,panel_name,ui, logger, row, col, drow, dcol):
        #super(image_panel, self).__init__()
        # enable some user interaction
        #fi.get_bindings.enable_all(True)

        self.logger = logger
        self.ui = ui
        self.dc = get_canvas_types()
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        #fi.set_callback('drag-drop', self.drop_file)
        #fi.set_callback('none-move',self.cursor_cb)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        fi.show_focus_indicator(True)
        self.fitsimage = fi
        
        bd = fi.get_bindings()
        bd.enable_all(True)
        self.cutout = fi.get_widget()
        self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.fitsimage = fi

        canvas = self.dc.DrawingCanvas()
        canvas.enable_draw(True)
        canvas.enable_edit(True)
        canvas.set_drawtype('rectangle', color='lightblue')
        canvas.set_surface(fi)
        #canvas.rectangle(.5,.5,10,0)
        self.canvas = canvas
        # add canvas to view
        #fi.add(canvas)
        private_canvas = fi.get_canvas()
        private_canvas.add(canvas)
        canvas.register_for_cursor_drawing(fi)
        #canvas.add_callback('draw-event', self.draw_cb)
        canvas.set_draw_mode('draw')
        canvas.ui_set_active(True)
        self.canvas = canvas

        self.drawtypes = canvas.get_drawtypes()
        self.drawtypes.sort()
    def load_image(self, imagearray):
        #self.fitsimage.set_image(imagearray)
        self.fitsimage.set_data(imagearray)

    def load_file(self, filepath):
        image = load_data(filepath, logger=self.logger)
        self.image_wcs = WCS(filepath)
        self.fitsimage.set_image(image)
        #self.setWindowTitle(filepath)


class mplWindow(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super(mplWindow, self).__init__(parent)



        # a figure instance to plot on
        self.figure = Figure(figsize=(2,2))
        self.axes = self.figure.add_subplot(111)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
        self.setParent(parent)
        
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        # set the layout
        #self.plot()
        #parent.addWidget(self.canvas)
        #self.show()
        if parent:
            parent.addWidget(self.canvas,1,2,4,1)
        else:
            layout = QtWidgets.QGridLayout(parent)
            layout.addWidget(self.canvas,1,2,4,1)
        self.setLayout(layout)
  
    def show(self):
        self.axes.plot(np.arange(10))
        self.canvas.draw()

