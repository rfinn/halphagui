#!/usr/bin/env python



# add a frame on top to show
# cutout, mask, model, residual



if __name__ == "__main__":
    #catalog = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    #from halphamain import cutout_image
    #from halphamain import cutout_image
    logger = log.get_logger("masklog", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    MainWindow = QtWidgets.QWidget()
    ui = galfitwindow(MainWindow, logger)
    #ui.setupUi(MainWindow)
    #ui.test()

    MainWindow.show()
    sys.exit(app.exec_())

    
