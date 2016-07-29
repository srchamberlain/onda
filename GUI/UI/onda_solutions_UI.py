# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'OndaCrystallographyGUI.ui'
#
# Created: Tue Dec 15 21:57:42 2015
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!
#
#
#    Added July 7, 2016
#    By Sarah Chamberlain, BioXFEL contribution

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_mainWindow(object):
    def setupUi(self, mainWindow):
        mainWindow.setObjectName(_fromUtf8("mainWindow"))
        mainWindow.resize(1200, 800)
        self.centralwidget = QtGui.QWidget(mainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout0 = QtGui.QVBoxLayout()
        self.verticalLayout0.setObjectName(_fromUtf8("verticalLayout0"))
        self.splitter0 = QtGui.QSplitter(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.splitter0.sizePolicy().hasHeightForWidth())
        self.splitter0.setSizePolicy(sizePolicy)
        self.splitter0.setOrientation(QtCore.Qt.Horizontal)
        self.splitter0.setObjectName(_fromUtf8("splitter0"))
        self.imageView = ImageView(self.splitter0)
        self.imageView.setObjectName(_fromUtf8("imageView"))
        self.splitter1 = QtGui.QSplitter(self.splitter0)
        self.splitter1.setOrientation(QtCore.Qt.Vertical)
        self.splitter1.setObjectName(_fromUtf8("splitter1"))
        self.splitter2 = QtGui.QSplitter(self.splitter0)    #Added
        self.splitter2.setOrientation(QtCore.Qt.Vertical)   #
        self.splitter2.setObjectName(_fromUtf8("splitter2"))    #
        self.unscaledPlotWidget = PlotWidget(self.splitter1)
        self.unscaledPlotWidget.setObjectName(_fromUtf8("unscaledPlotWidget"))
        self.scaledPlotWidget = PlotWidget(self.splitter1)
        self.scaledPlotWidget.setObjectName(_fromUtf8("qscaledPlotWidget"))
        self.sumPlotWidget = PlotWidget(self.splitter1)     #
        self.sumPlotWidget.setObjectName(_fromUtf8(("sumPlotWidget")))       #
        ## Not currently expanded to include specilized SAXS calculations
        self.RgPlotWidget = PlotWidget(self.splitter2) #Added
        self.RgPlotWidget.setObjectName(_fromUtf8("RgPlotWidget")) #
        self.intensityPlotWidget = PlotWidget(self.splitter2)
        self.intensityPlotWidget.setObjectName(_fromUtf8("intensityPlotWidget"))
        self.verticalLayout0.addWidget(self.splitter0)
        self.horizontalLayout0 = QtGui.QHBoxLayout()
        self.horizontalLayout0.setObjectName(_fromUtf8("horizontalLayout0"))
        self.savePlotButton = QtGui.QPushButton(self.centralwidget)  #Button for saving Cumulative average plot
        self.savePlotButton.setObjectName(_fromUtf8("savePlotButton"))
        self.horizontalLayout0.addWidget(self.savePlotButton)
        self.resetPlotsButton = QtGui.QPushButton(self.centralwidget) #Button for reseting Plots
        self.resetPlotsButton.setObjectName(_fromUtf8("resetPlotsButton"))
        self.horizontalLayout0.addWidget(self.resetPlotsButton)
        self.comparePlotsButton = QtGui.QPushButton(self.centralwidget) #Button to toggle on/off compare plot
        self.comparePlotsButton.setObjectName(_fromUtf8("comparePlotsButton"))
        self.horizontalLayout0.addWidget(self.comparePlotsButton)
        self.xAxisTogglePlotsButton = QtGui.QPushButton(self.centralwidget) #Button to toggle xAxis
        self.xAxisTogglePlotsButton.setObjectName(_fromUtf8("xAxisTogglePlotsButton"))
        self.horizontalLayout0.addWidget(self.xAxisTogglePlotsButton)
        self.NumPlotLabel = QtGui.QLabel(self.centralwidget) # Number of Profiles in cumulative average profile
        self.NumPlotLabel.setObjectName(_fromUtf8("NumPlotLabel"))
        self.horizontalLayout0.addWidget(self.NumPlotLabel)
        self.NumTotalLabel = QtGui.QLabel(self.centralwidget) # Num of Profiles tested w/ std. dev for cumulative average profile
        self.NumTotalLabel.setObjectName(_fromUtf8("NumTotalLabel"))
        self.horizontalLayout0.addWidget(self.NumTotalLabel)
        self.percentPlotLabel = QtGui.QLabel(self.centralwidget) # Percent of Profiles Plotted in Cum. plot
        self.percentPlotLabel.setObjectName(_fromUtf8("percentPlotLabel"))
        self.horizontalLayout0.addWidget(self.percentPlotLabel)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout0.addItem(spacerItem)
        self.verticalLayout0.addLayout(self.horizontalLayout0)
        self.gridLayout.addLayout(self.verticalLayout0, 0, 0, 1, 1)
        mainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(mainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        mainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(mainWindow)
        QtCore.QMetaObject.connectSlotsByName(mainWindow)

    def retranslateUi(self, mainWindow):
        mainWindow.setWindowTitle(_translate("mainWindow", "OnDA Solutions", None))
        # Defines button labels in GUI
        self.savePlotButton.setText(_translate("mainWindow", "Save Cumulative Plot", None))
        self.resetPlotsButton.setText(_translate("mainWindow", "Reset Plots", None))
        self.comparePlotsButton.setText(_translate("mainWindow", "Compare Plots", None))
        self.xAxisTogglePlotsButton.setText(_translate("mainWindow", "x-axis Toggle", None))
        self.NumPlotLabel.setText(_translate("mainWindow", "Number Averaged: -", None))
        self.NumTotalLabel.setText(_translate("mainWindow", "Number Processed: -", None))
        self.percentPlotLabel.setText(_translate("mainWindow", "% Plotted: -", None))

from pyqtgraph import ImageView, PlotWidget
