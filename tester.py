from PyQt5.QtCore import QDate, QTime, QDateTime, Qt

now = QDate.currentDate()

print(now.toString(Qt.ISODate))

time = QTime.currentTime()

print(time.toString(Qt.DefaultLocaleShortDate))
name = "Alex"
surname = "Zaplatil"
name_of_sim = "Horiton2"
f = open(name_of_sim + now.toString(Qt.ISODate) + "-" + time.toString(Qt.DefaultLocaleShortDate)[:2] +
                 "-" + time.toString(Qt.DefaultLocaleShortDate)[-2:] + ".txt", "w+")
