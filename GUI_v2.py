# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1080, 754)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.menuWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.menuWidget.setEnabled(True)
        self.menuWidget.setGeometry(QtCore.QRect(0, 0, 1080, 720))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.menuWidget.setFont(font)
        self.menuWidget.setObjectName("menuWidget")
        self.Simulation = QtWidgets.QWidget()
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Simulation.setFont(font)
        self.Simulation.setObjectName("Simulation")
        self.track_photo_sim = QtWidgets.QLabel(self.Simulation)
        self.track_photo_sim.setGeometry(QtCore.QRect(620, 80, 431, 331))
        self.track_photo_sim.setText("")
        self.track_photo_sim.setPixmap(QtGui.QPixmap("Photos/Capture.PNG"))
        self.track_photo_sim.setScaledContents(True)
        self.track_photo_sim.setObjectName("track_photo_sim")
        self.infoBox = QtWidgets.QGroupBox(self.Simulation)
        self.infoBox.setGeometry(QtCore.QRect(50, 80, 391, 341))
        self.infoBox.setObjectName("infoBox")
        self.author = QtWidgets.QLabel(self.infoBox)
        self.author.setGeometry(QtCore.QRect(10, 30, 121, 21))
        self.author.setObjectName("author")
        self.date = QtWidgets.QLabel(self.infoBox)
        self.date.setGeometry(QtCore.QRect(10, 70, 71, 31))
        self.date.setObjectName("date")
        self.name_of_sim = QtWidgets.QLabel(self.infoBox)
        self.name_of_sim.setGeometry(QtCore.QRect(10, 120, 181, 16))
        self.name_of_sim.setObjectName("name_of_sim")
        self.note = QtWidgets.QLabel(self.infoBox)
        self.note.setGeometry(QtCore.QRect(10, 160, 55, 16))
        self.note.setObjectName("note")
        self.dateTimeEdit = QtWidgets.QDateTimeEdit(self.infoBox)
        self.dateTimeEdit.setGeometry(QtCore.QRect(140, 70, 194, 22))
        self.dateTimeEdit.setObjectName("dateTimeEdit")
        self.name_lineEdit = QtWidgets.QLineEdit(self.infoBox)
        self.name_lineEdit.setGeometry(QtCore.QRect(140, 30, 221, 22))
        self.name_lineEdit.setObjectName("name_lineEdit")
        self.name_of_sim_Edit = QtWidgets.QLineEdit(self.infoBox)
        self.name_of_sim_Edit.setGeometry(QtCore.QRect(200, 120, 161, 22))
        self.name_of_sim_Edit.setObjectName("name_of_sim_Edit")
        self.write_box = QtWidgets.QTextEdit(self.infoBox)
        self.write_box.setGeometry(QtCore.QRect(10, 180, 361, 151))
        self.write_box.setObjectName("write_box")
        self.simulateButton = QtWidgets.QPushButton(self.Simulation)
        self.simulateButton.setGeometry(QtCore.QRect(50, 470, 191, 81))
        self.simulateButton.setObjectName("simulateButton")
        self.data_download_Box = QtWidgets.QComboBox(self.Simulation)
        self.data_download_Box.setGeometry(QtCore.QRect(340, 470, 201, 81))
        self.data_download_Box.setObjectName("data_download_Box")
        self.data_download_Box.addItem("")
        self.data_download_Box.addItem("")
        self.exportButton = QtWidgets.QPushButton(self.Simulation)
        self.exportButton.setGeometry(QtCore.QRect(550, 470, 171, 81))
        self.exportButton.setObjectName("exportButton")
        self.menuWidget.addTab(self.Simulation, "")
        self.Track = QtWidgets.QWidget()
        self.Track.setObjectName("Track")
        self.loadtrackButton = QtWidgets.QPushButton(self.Track)
        self.loadtrackButton.setGeometry(QtCore.QRect(50, 480, 171, 51))
        self.loadtrackButton.setObjectName("loadtrackButton")
        self.select_track = QtWidgets.QLabel(self.Track)
        self.select_track.setGeometry(QtCore.QRect(50, 80, 151, 31))
        self.select_track.setObjectName("select_track")
        self.track_photo = QtWidgets.QLabel(self.Track)
        self.track_photo.setGeometry(QtCore.QRect(500, 130, 391, 331))
        self.track_photo.setText("")
        self.track_photo.setPixmap(QtGui.QPixmap("Photos/Capture.PNG"))
        self.track_photo.setScaledContents(True)
        self.track_photo.setObjectName("track_photo")
        self.scrollArea = QtWidgets.QScrollArea(self.Track)
        self.scrollArea.setGeometry(QtCore.QRect(50, 120, 301, 331))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.scrollArea.setFont(font)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 299, 329))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.scrollAreaWidgetContents.setFont(font)
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.listWidget = QtWidgets.QListWidget(self.scrollAreaWidgetContents)
        self.listWidget.setGeometry(QtCore.QRect(0, 0, 311, 331))
        self.listWidget.setObjectName("listWidget")
        item = QtWidgets.QListWidgetItem()
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listWidget.addItem(item)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.menuWidget.addTab(self.Track, "")
        self.Basic = QtWidgets.QWidget()
        self.Basic.setObjectName("Basic")
        self.Basic_frame = QtWidgets.QFrame(self.Basic)
        self.Basic_frame.setGeometry(QtCore.QRect(0, 0, 1041, 581))
        self.Basic_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Basic_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.Basic_frame.setObjectName("Basic_frame")
        self.basic_tableWidget = QtWidgets.QTableWidget(self.Basic_frame)
        self.basic_tableWidget.setGeometry(QtCore.QRect(10, 10, 370, 551))
        self.basic_tableWidget.setObjectName("basic_tableWidget")
        self.basic_tableWidget.setColumnCount(2)
        self.basic_tableWidget.setRowCount(6)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setHorizontalHeaderItem(1, item)
        self.basic_textEdit = QtWidgets.QTextEdit(self.Basic_frame)
        self.basic_textEdit.setGeometry(QtCore.QRect(530, 60, 451, 231))
        self.basic_textEdit.setObjectName("basic_textEdit")
        self.explanation_basic = QtWidgets.QLabel(self.Basic_frame)
        self.explanation_basic.setGeometry(QtCore.QRect(530, 20, 131, 31))
        self.explanation_basic.setObjectName("explanation_basic")
        self.save_basic_Button = QtWidgets.QPushButton(self.Basic_frame)
        self.save_basic_Button.setGeometry(QtCore.QRect(530, 320, 171, 81))
        self.save_basic_Button.setObjectName("save_basic_Button")
        self.load_basic_Button = QtWidgets.QPushButton(self.Basic_frame)
        self.load_basic_Button.setGeometry(QtCore.QRect(720, 320, 191, 81))
        self.load_basic_Button.setObjectName("load_basic_Button")
        self.menuWidget.addTab(self.Basic, "")
        self.Suspension = QtWidgets.QWidget()
        self.Suspension.setObjectName("Suspension")
        self.suspension_frame = QtWidgets.QFrame(self.Suspension)
        self.suspension_frame.setGeometry(QtCore.QRect(0, 0, 1041, 581))
        self.suspension_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.suspension_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.suspension_frame.setObjectName("suspension_frame")
        self.suspension_tableWidget = QtWidgets.QTableWidget(self.suspension_frame)
        self.suspension_tableWidget.setGeometry(QtCore.QRect(10, 10, 370, 551))
        self.suspension_tableWidget.setObjectName("suspension_tableWidget")
        self.suspension_tableWidget.setColumnCount(2)
        self.suspension_tableWidget.setRowCount(1)
        item = QtWidgets.QTableWidgetItem()
        self.suspension_tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.suspension_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.suspension_tableWidget.setHorizontalHeaderItem(1, item)
        self.suspension_textEdit = QtWidgets.QTextEdit(self.suspension_frame)
        self.suspension_textEdit.setGeometry(QtCore.QRect(530, 60, 451, 231))
        self.suspension_textEdit.setObjectName("suspension_textEdit")
        self.explanation_basic_2 = QtWidgets.QLabel(self.suspension_frame)
        self.explanation_basic_2.setGeometry(QtCore.QRect(530, 20, 131, 31))
        self.explanation_basic_2.setObjectName("explanation_basic_2")
        self.save_suspension_Button = QtWidgets.QPushButton(self.suspension_frame)
        self.save_suspension_Button.setGeometry(QtCore.QRect(530, 320, 171, 81))
        self.save_suspension_Button.setObjectName("save_suspension_Button")
        self.load_suspension_Button = QtWidgets.QPushButton(self.suspension_frame)
        self.load_suspension_Button.setGeometry(QtCore.QRect(720, 320, 191, 81))
        self.load_suspension_Button.setObjectName("load_suspension_Button")
        self.menuWidget.addTab(self.Suspension, "")
        self.Aerodynamics = QtWidgets.QWidget()
        self.Aerodynamics.setObjectName("Aerodynamics")
        self.aerodynamics_frame = QtWidgets.QFrame(self.Aerodynamics)
        self.aerodynamics_frame.setGeometry(QtCore.QRect(0, 0, 1041, 581))
        self.aerodynamics_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.aerodynamics_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.aerodynamics_frame.setObjectName("aerodynamics_frame")
        self.aerodynamics_tableWidget = QtWidgets.QTableWidget(self.aerodynamics_frame)
        self.aerodynamics_tableWidget.setGeometry(QtCore.QRect(10, 10, 370, 551))
        self.aerodynamics_tableWidget.setObjectName("aerodynamics_tableWidget")
        self.aerodynamics_tableWidget.setColumnCount(2)
        self.aerodynamics_tableWidget.setRowCount(6)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setHorizontalHeaderItem(1, item)
        self.aerodynamics_textEdit = QtWidgets.QTextEdit(self.aerodynamics_frame)
        self.aerodynamics_textEdit.setGeometry(QtCore.QRect(530, 60, 451, 231))
        self.aerodynamics_textEdit.setObjectName("aerodynamics_textEdit")
        self.explanation_aerodnyamics = QtWidgets.QLabel(self.aerodynamics_frame)
        self.explanation_aerodnyamics.setGeometry(QtCore.QRect(530, 20, 131, 31))
        self.explanation_aerodnyamics.setObjectName("explanation_aerodnyamics")
        self.save_aerodynamics_Button = QtWidgets.QPushButton(self.aerodynamics_frame)
        self.save_aerodynamics_Button.setGeometry(QtCore.QRect(530, 320, 171, 81))
        self.save_aerodynamics_Button.setObjectName("save_aerodynamics_Button")
        self.load_aerodynamics_Button = QtWidgets.QPushButton(self.aerodynamics_frame)
        self.load_aerodynamics_Button.setGeometry(QtCore.QRect(720, 320, 191, 81))
        self.load_aerodynamics_Button.setObjectName("load_aerodynamics_Button")
        self.menuWidget.addTab(self.Aerodynamics, "")
        self.Drivetrain = QtWidgets.QWidget()
        self.Drivetrain.setObjectName("Drivetrain")
        self.drivetrain_frame = QtWidgets.QFrame(self.Drivetrain)
        self.drivetrain_frame.setGeometry(QtCore.QRect(0, 0, 1041, 581))
        self.drivetrain_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.drivetrain_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.drivetrain_frame.setObjectName("drivetrain_frame")
        self.drivetrain_tableWidget = QtWidgets.QTableWidget(self.drivetrain_frame)
        self.drivetrain_tableWidget.setGeometry(QtCore.QRect(10, 10, 370, 551))
        self.drivetrain_tableWidget.setObjectName("drivetrain_tableWidget")
        self.drivetrain_tableWidget.setColumnCount(2)
        self.drivetrain_tableWidget.setRowCount(5)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setHorizontalHeaderItem(1, item)
        self.drivetrain_textEdit = QtWidgets.QTextEdit(self.drivetrain_frame)
        self.drivetrain_textEdit.setGeometry(QtCore.QRect(530, 60, 451, 231))
        self.drivetrain_textEdit.setObjectName("drivetrain_textEdit")
        self.explanation_drivetrain = QtWidgets.QLabel(self.drivetrain_frame)
        self.explanation_drivetrain.setGeometry(QtCore.QRect(530, 20, 131, 31))
        self.explanation_drivetrain.setObjectName("explanation_drivetrain")
        self.save_drivetrain_Button = QtWidgets.QPushButton(self.drivetrain_frame)
        self.save_drivetrain_Button.setGeometry(QtCore.QRect(530, 320, 171, 81))
        self.save_drivetrain_Button.setObjectName("save_drivetrain_Button")
        self.load_drivetrain_Button = QtWidgets.QPushButton(self.drivetrain_frame)
        self.load_drivetrain_Button.setGeometry(QtCore.QRect(720, 320, 191, 81))
        self.load_drivetrain_Button.setObjectName("load_drivetrain_Button")
        self.menuWidget.addTab(self.Drivetrain, "")
        self.Electronics = QtWidgets.QWidget()
        self.Electronics.setObjectName("Electronics")
        self.electronics_frame = QtWidgets.QFrame(self.Electronics)
        self.electronics_frame.setGeometry(QtCore.QRect(0, 0, 1041, 581))
        self.electronics_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.electronics_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.electronics_frame.setObjectName("electronics_frame")
        self.electronics_tableWidget = QtWidgets.QTableWidget(self.electronics_frame)
        self.electronics_tableWidget.setGeometry(QtCore.QRect(10, 10, 370, 551))
        self.electronics_tableWidget.setObjectName("electronics_tableWidget")
        self.electronics_tableWidget.setColumnCount(2)
        self.electronics_tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.electronics_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.electronics_tableWidget.setHorizontalHeaderItem(1, item)
        self.electronics_textEdit = QtWidgets.QTextEdit(self.electronics_frame)
        self.electronics_textEdit.setGeometry(QtCore.QRect(530, 60, 451, 231))
        self.electronics_textEdit.setObjectName("electronics_textEdit")
        self.explanation_electronics = QtWidgets.QLabel(self.electronics_frame)
        self.explanation_electronics.setGeometry(QtCore.QRect(530, 20, 131, 31))
        self.explanation_electronics.setObjectName("explanation_electronics")
        self.save_electronics_Button = QtWidgets.QPushButton(self.electronics_frame)
        self.save_electronics_Button.setGeometry(QtCore.QRect(530, 320, 171, 81))
        self.save_electronics_Button.setObjectName("save_electronics_Button")
        self.load_electronics_Button = QtWidgets.QPushButton(self.electronics_frame)
        self.load_electronics_Button.setGeometry(QtCore.QRect(720, 320, 191, 81))
        self.load_electronics_Button.setObjectName("load_electronics_Button")
        self.menuWidget.addTab(self.Electronics, "")
        self.View = QtWidgets.QWidget()
        self.View.setObjectName("View")
        self.menuWidget.addTab(self.View, "")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1080, 26))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionNew = QtWidgets.QAction(MainWindow)
        self.actionNew.setObjectName("actionNew")
        self.actionSave = QtWidgets.QAction(MainWindow)
        self.actionSave.setObjectName("actionSave")
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionExit = QtWidgets.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.menuFile.addAction(self.actionNew)
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        self.menuWidget.setCurrentIndex(6)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.infoBox.setTitle(_translate("MainWindow", "INFO"))
        self.author.setText(_translate("MainWindow", "Author"))
        self.date.setText(_translate("MainWindow", "Date"))
        self.name_of_sim.setText(_translate("MainWindow", "Name of Simulation"))
        self.note.setText(_translate("MainWindow", "Notes"))
        self.name_lineEdit.setText(_translate("MainWindow", "Name and Surname"))
        self.name_of_sim_Edit.setText(_translate("MainWindow", "Name"))
        self.write_box.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:12pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-style:italic;\">write changes, notes</span></p></body></html>"))
        self.simulateButton.setText(_translate("MainWindow", "Simulate"))
        self.data_download_Box.setItemText(0, _translate("MainWindow", "CSV- basic"))
        self.data_download_Box.setItemText(1, _translate("MainWindow", "CSV- advanced"))
        self.exportButton.setText(_translate("MainWindow", "Export"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Simulation), _translate("MainWindow", "Simulation"))
        self.loadtrackButton.setText(_translate("MainWindow", "Load track"))
        self.select_track.setText(_translate("MainWindow", "Select track"))
        __sortingEnabled = self.listWidget.isSortingEnabled()
        self.listWidget.setSortingEnabled(False)
        item = self.listWidget.item(0)
        item.setText(_translate("MainWindow", "Acceleration"))
        item = self.listWidget.item(1)
        item.setText(_translate("MainWindow", "Skidpad"))
        item = self.listWidget.item(2)
        item.setText(_translate("MainWindow", "FSG 2012 Autocross"))
        self.listWidget.setSortingEnabled(__sortingEnabled)
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Track), _translate("MainWindow", "Track"))
        item = self.basic_tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "Mass"))
        item = self.basic_tableWidget.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "CG in y"))
        item = self.basic_tableWidget.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "CG in z"))
        item = self.basic_tableWidget.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "Track widht"))
        item = self.basic_tableWidget.verticalHeaderItem(4)
        item.setText(_translate("MainWindow", "Wheelbase"))
        item = self.basic_tableWidget.verticalHeaderItem(5)
        item.setText(_translate("MainWindow", "W"))
        item = self.basic_tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Value"))
        item = self.basic_tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Unit"))
        self.explanation_basic.setText(_translate("MainWindow", "Explanation"))
        self.save_basic_Button.setText(_translate("MainWindow", "Save"))
        self.load_basic_Button.setText(_translate("MainWindow", "Load"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Basic), _translate("MainWindow", "Basic"))
        item = self.suspension_tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "Coefficient of Friction"))
        item = self.suspension_tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Value"))
        item = self.suspension_tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Unit"))
        self.explanation_basic_2.setText(_translate("MainWindow", "Explanation"))
        self.save_suspension_Button.setText(_translate("MainWindow", "Save"))
        self.load_suspension_Button.setText(_translate("MainWindow", "Load"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Suspension), _translate("MainWindow", "Suspension"))
        item = self.aerodynamics_tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "Air density"))
        item = self.aerodynamics_tableWidget.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "Frontal area"))
        item = self.aerodynamics_tableWidget.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "Coefficient of DF"))
        item = self.aerodynamics_tableWidget.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "Coefficient of Drag"))
        item = self.aerodynamics_tableWidget.verticalHeaderItem(4)
        item.setText(_translate("MainWindow", "CoP in y"))
        item = self.aerodynamics_tableWidget.verticalHeaderItem(5)
        item.setText(_translate("MainWindow", "CoP in z"))
        item = self.aerodynamics_tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Value"))
        item = self.aerodynamics_tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Unit"))
        self.explanation_aerodnyamics.setText(_translate("MainWindow", "Explanation"))
        self.save_aerodynamics_Button.setText(_translate("MainWindow", "Save"))
        self.load_aerodynamics_Button.setText(_translate("MainWindow", "Load"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Aerodynamics), _translate("MainWindow", "Aerodynamics"))
        item = self.drivetrain_tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "Tire radius"))
        item = self.drivetrain_tableWidget.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "Max Power"))
        item = self.drivetrain_tableWidget.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "Max RPM"))
        item = self.drivetrain_tableWidget.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "Gear ratio"))
        item = self.drivetrain_tableWidget.verticalHeaderItem(4)
        item.setText(_translate("MainWindow", "Max torque"))
        item = self.drivetrain_tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Value"))
        item = self.drivetrain_tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Unit"))
        self.explanation_drivetrain.setText(_translate("MainWindow", "Explanation"))
        self.save_drivetrain_Button.setText(_translate("MainWindow", "Save"))
        self.load_drivetrain_Button.setText(_translate("MainWindow", "Load"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Drivetrain), _translate("MainWindow", "Drivetrain"))
        item = self.electronics_tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Value"))
        item = self.electronics_tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Unit"))
        self.explanation_electronics.setText(_translate("MainWindow", "Explanation"))
        self.save_electronics_Button.setText(_translate("MainWindow", "Save"))
        self.load_electronics_Button.setText(_translate("MainWindow", "Load"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.Electronics), _translate("MainWindow", "Electronics"))
        self.menuWidget.setTabText(self.menuWidget.indexOf(self.View), _translate("MainWindow", "View"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionNew.setText(_translate("MainWindow", "New"))
        self.actionSave.setText(_translate("MainWindow", "Save"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
