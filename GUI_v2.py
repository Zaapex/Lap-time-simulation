# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox, QTableWidgetItem, QFileDialog
from PyQt5.QtCore import QDate, QTime, Qt
import pandas as pd
from tester import *
from Code import *

selected_track = ["None"]
selected_settings = ["Svarog data"]
df = pd.read_csv(selected_settings[0])

# time
now = QDate.currentDate()
time = QTime.currentTime()

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

        # simulation widget
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
        #self.dateTimeEdit = QtWidgets.QDateTimeEdit(self.infoBox)
        self.dateTimeEdit = QtWidgets.QLabel(self.infoBox)
        self.dateTimeEdit.setGeometry(QtCore.QRect(140, 70, 194, 22))
        self.dateTimeEdit.setObjectName("dateTimeEdit")
        self.name_lineEdit = QtWidgets.QLineEdit(self.infoBox)
        self.name_lineEdit.setGeometry(QtCore.QRect(140, 30, 221, 22))
        self.name_lineEdit.setObjectName("name_lineEdit")
        self.name_of_sim_Edit = QtWidgets.QLineEdit(self.infoBox)
        self.name_of_sim_Edit.setGeometry(QtCore.QRect(200, 120, 161, 22))
        self.name_of_sim_Edit.setObjectName("name_of_sim_Edit")
        self.write_box = QtWidgets.QLineEdit(self.infoBox)
        self.write_box.setGeometry(QtCore.QRect(10, 180, 361, 151))
        self.write_box.setObjectName("write_box")
        self.simulateButton = QtWidgets.QPushButton(self.Simulation)
        self.simulateButton.setGeometry(QtCore.QRect(50, 470, 191, 81))
        self.simulateButton.setObjectName("simulateButton")
        self.simulateButton.clicked.connect(self.clickedinfo)  # ADDED
        self.data_download_Box = QtWidgets.QComboBox(self.Simulation)
        self.data_download_Box.setGeometry(QtCore.QRect(340, 470, 201, 81))
        self.data_download_Box.setObjectName("data_download_Box")
        self.data_download_Box.addItem("")
        self.data_download_Box.addItem("")
        self.exportButton = QtWidgets.QPushButton(self.Simulation)
        self.exportButton.setGeometry(QtCore.QRect(550, 470, 171, 81))
        self.exportButton.setObjectName("exportButton")
        self.menuWidget.addTab(self.Simulation, "")

        # track widget
        self.Track = QtWidgets.QWidget()
        self.Track.setObjectName("Track")
        self.loadtrackButton = QtWidgets.QPushButton(self.Track)
        self.loadtrackButton.setGeometry(QtCore.QRect(50, 480, 171, 51))
        self.loadtrackButton.setObjectName("loadtrackButton")
        self.loadtrackButton.clicked.connect(self.load_track_clicked)
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
        self.listWidget.itemClicked.connect(self.track_cliked)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.menuWidget.addTab(self.Track, "")

        # basic widget
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
        mass_index = df.loc[df['Parameter'] == "Mass"].index[0]  # get mass value
        self.basic_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][mass_index])))  # write in the table
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(1, item)
        cg_y_index = df.loc[df['Parameter'] == "CG in y"].index[0]
        self.basic_tableWidget.setItem(1, 0, QTableWidgetItem(str(df["Value"][cg_y_index])))
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(2, item)
        cg_z_index = df.loc[df['Parameter'] == "CG in z"].index[0]
        self.basic_tableWidget.setItem(2, 0, QTableWidgetItem(str(df["Value"][cg_z_index])))
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(3, item)
        track_width_index = df.loc[df['Parameter'] == "Track width"].index[0]
        self.basic_tableWidget.setItem(3, 0, QTableWidgetItem(str(df["Value"][track_width_index])))
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(4, item)
        wheelbase_index = df.loc[df['Parameter'] == "Wheelbase"].index[0]
        self.basic_tableWidget.setItem(4, 0, QTableWidgetItem(str(df["Value"][wheelbase_index])))
        item = QtWidgets.QTableWidgetItem()
        self.basic_tableWidget.setVerticalHeaderItem(5, item)
        w_index = df.loc[df['Parameter'] == "W"].index[0]
        self.basic_tableWidget.setItem(5, 0, QTableWidgetItem(str(df["Value"][w_index])))
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
        self.load_basic_Button.clicked.connect(self.get_formula_data)
        self.menuWidget.addTab(self.Basic, "")

        # suspension widget
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
        coef_friction_index = df.loc[df['Parameter'] == "Coefficient of Friction"].index[0]
        self.suspension_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][coef_friction_index])))
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
        self.load_suspension_Button.clicked.connect(self.get_formula_data)
        self.menuWidget.addTab(self.Suspension, "")

        # aerodynamics widget
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
        air_index = df.loc[df['Parameter'] == "Air density"].index[0]
        self.aerodynamics_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][air_index])))
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(1, item)
        f_area_index = df.loc[df['Parameter'] == "Frontal area"].index[0]
        self.aerodynamics_tableWidget.setItem(1, 0, QTableWidgetItem(str(df["Value"][f_area_index])))
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(2, item)
        coef_df_index = df.loc[df['Parameter'] == "Coefficient of DF"].index[0]
        self.aerodynamics_tableWidget.setItem(2, 0, QTableWidgetItem(str(df["Value"][coef_df_index])))
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(3, item)
        coef_drag_index = df.loc[df['Parameter'] == "Coefficient of Drag"].index[0]
        self.aerodynamics_tableWidget.setItem(3, 0, QTableWidgetItem(str(df["Value"][coef_drag_index])))
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(4, item)
        cop_y_index = df.loc[df['Parameter'] == "CoP in y"].index[0]
        self.aerodynamics_tableWidget.setItem(4, 0, QTableWidgetItem(str(df["Value"][cop_y_index])))
        item = QtWidgets.QTableWidgetItem()
        self.aerodynamics_tableWidget.setVerticalHeaderItem(5, item)
        cop_z_index = df.loc[df['Parameter'] == "CoP in z"].index[0]
        self.aerodynamics_tableWidget.setItem(5, 0, QTableWidgetItem(str(df["Value"][cop_z_index])))
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
        self.load_aerodynamics_Button.clicked.connect(self.get_formula_data)
        self.menuWidget.addTab(self.Aerodynamics, "")

        # Drivetrain widget
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
        tire_radius_index = df.loc[df['Parameter'] == "Tire radius"].index[0]
        self.drivetrain_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][tire_radius_index])))
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(1, item)
        max_power_index = df.loc[df['Parameter'] == "Max Power"].index[0]
        self.drivetrain_tableWidget.setItem(1, 0, QTableWidgetItem(str(df["Value"][max_power_index])))
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(2, item)
        max_rpm_index = df.loc[df['Parameter'] == "Max RPM"].index[0]
        self.drivetrain_tableWidget.setItem(2, 0, QTableWidgetItem(str(df["Value"][max_rpm_index])))
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(3, item)
        gear_ratio_index = df.loc[df['Parameter'] == "Gear ratio"].index[0]
        self.drivetrain_tableWidget.setItem(3, 0, QTableWidgetItem(str(df["Value"][gear_ratio_index])))
        item = QtWidgets.QTableWidgetItem()
        self.drivetrain_tableWidget.setVerticalHeaderItem(4, item)
        max_torque_index = df.loc[df['Parameter'] == "Max torque"].index[0]
        self.drivetrain_tableWidget.setItem(4, 0, QTableWidgetItem(str(df["Value"][max_torque_index])))
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
        self.load_drivetrain_Button.clicked.connect(self.get_formula_data)
        self.menuWidget.addTab(self.Drivetrain, "")

        # Electronics widget
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
        self.load_electronics_Button.clicked.connect(self.get_formula_data)
        self.menuWidget.addTab(self.Electronics, "")

        # view widget
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
        self.dateTimeEdit.setText(str(now.toString(Qt.ISODate)) + " " + time.toString(Qt.DefaultLocaleShortDate))
        self.name_of_sim.setText(_translate("MainWindow", "Name of Simulation"))
        self.note.setText(_translate("MainWindow", "Notes"))
        self.name_lineEdit.setText(_translate("MainWindow", "Name and Surname"))
        self.name_of_sim_Edit.setText(_translate("MainWindow", "Name"))
        self.write_box.setText(_translate("MainWindow", "Write notes"))
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

    def clickedinfo(self):
        author = self.name_lineEdit.text()
        name_of_sim = self.name_of_sim_Edit.text()
        notes = self.note.text()

        f = open(name_of_sim + now.toString(Qt.ISODate) + "-" + time.toString(Qt.DefaultLocaleShortDate)[:2] +
                 "-" + time.toString(Qt.DefaultLocaleShortDate)[-2:] + ".txt", "w+")

        f.write("Author: " + str(author) + "\n")
        f.write("Date: " + str(now.toString(Qt.ISODate)) + "\n")
        f.write("Name of Simulation: " + str(name_of_sim) + "\n")
        f.write("Notes: " + str(notes) + "\n")

        lap_time = lap_time_simulation(selected_track[0], selected_settings[0])

        for line in range(len(df.index)):
            f.write(df["Parameter"][line] + ": " + df["Value"][line] + "\n")

        f.write("Total time: " + str(lap_time) + "\n")

        msg = QMessageBox()
        msg.setText("Simulation is finished")
        msg.setWindowTitle("Simulation")

        x = msg.exec_()

    def track_cliked(self, item):
        track = item.text()
        self.track_photo.setPixmap(QtGui.QPixmap("Photos/" + str(track) + ".jpg"))
        selected_track[0] = track

    def load_track_clicked(self):
        msg = QMessageBox()
        msg.setText("You selected track: " + selected_track[0])
        msg.setWindowTitle("Track loaded")

        x = msg.exec_()

    def get_formula_data(self):
        fname = QFileDialog.getOpenFileName()
        selected_settings[0] = str(fname[0])
        # data changes when you click choose file

        df = pd.read_csv(selected_settings[0])
        mass_index = df.loc[df['Parameter'] == "Mass"].index[0]  # get mass value
        self.basic_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][mass_index])))  # write in the table
        cg_y_index = df.loc[df['Parameter'] == "CG in y"].index[0]
        self.basic_tableWidget.setItem(1, 0, QTableWidgetItem(str(df["Value"][cg_y_index])))
        cg_z_index = df.loc[df['Parameter'] == "CG in z"].index[0]
        self.basic_tableWidget.setItem(2, 0, QTableWidgetItem(str(df["Value"][cg_z_index])))
        track_width_index = df.loc[df['Parameter'] == "Track width"].index[0]
        self.basic_tableWidget.setItem(3, 0, QTableWidgetItem(str(df["Value"][track_width_index])))
        wheelbase_index = df.loc[df['Parameter'] == "Wheelbase"].index[0]
        self.basic_tableWidget.setItem(4, 0, QTableWidgetItem(str(df["Value"][wheelbase_index])))
        w_index = df.loc[df['Parameter'] == "W"].index[0]
        self.basic_tableWidget.setItem(5, 0, QTableWidgetItem(str(df["Value"][w_index])))

        coef_friction_index = df.loc[df['Parameter'] == "Coefficient of Friction"].index[0]
        self.suspension_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][coef_friction_index])))

        air_index = df.loc[df['Parameter'] == "Air density"].index[0]
        self.aerodynamics_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][air_index])))
        f_area_index = df.loc[df['Parameter'] == "Frontal area"].index[0]
        self.aerodynamics_tableWidget.setItem(1, 0, QTableWidgetItem(str(df["Value"][f_area_index])))
        coef_df_index = df.loc[df['Parameter'] == "Coefficient of DF"].index[0]
        self.aerodynamics_tableWidget.setItem(2, 0, QTableWidgetItem(str(df["Value"][coef_df_index])))
        coef_drag_index = df.loc[df['Parameter'] == "Coefficient of Drag"].index[0]
        self.aerodynamics_tableWidget.setItem(3, 0, QTableWidgetItem(str(df["Value"][coef_drag_index])))
        cop_y_index = df.loc[df['Parameter'] == "CoP in y"].index[0]
        self.aerodynamics_tableWidget.setItem(4, 0, QTableWidgetItem(str(df["Value"][cop_y_index])))
        cop_z_index = df.loc[df['Parameter'] == "CoP in z"].index[0]
        self.aerodynamics_tableWidget.setItem(5, 0, QTableWidgetItem(str(df["Value"][cop_z_index])))

        tire_radius_index = df.loc[df['Parameter'] == "Tire radius"].index[0]
        self.drivetrain_tableWidget.setItem(0, 0, QTableWidgetItem(str(df["Value"][tire_radius_index])))
        max_power_index = df.loc[df['Parameter'] == "Max Power"].index[0]
        self.drivetrain_tableWidget.setItem(1, 0, QTableWidgetItem(str(df["Value"][max_power_index])))
        max_rpm_index = df.loc[df['Parameter'] == "Max RPM"].index[0]
        self.drivetrain_tableWidget.setItem(2, 0, QTableWidgetItem(str(df["Value"][max_rpm_index])))
        gear_ratio_index = df.loc[df['Parameter'] == "Gear ratio"].index[0]
        self.drivetrain_tableWidget.setItem(3, 0, QTableWidgetItem(str(df["Value"][gear_ratio_index])))
        max_torque_index = df.loc[df['Parameter'] == "Max torque"].index[0]
        self.drivetrain_tableWidget.setItem(4, 0, QTableWidgetItem(str(df["Value"][max_torque_index])))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
