#ifndef FORM1_H
#define FORM1_H

#include <qvariant.h>


#include <Qt3Support/Q3MainWindow>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QPushButton>
#include <QtGui/QTabWidget>
#include <QtGui/QWidget>
#include <Qt3Support/Q3MimeSourceFactory>

class Ui_meshGenQt
{
public:
    QWidget *widget;
    QPushButton *saveButton;
    QPushButton *quitButton;
    QTabWidget *tabWidget;
    QWidget *tab;
    QWidget *tab1;
    QWidget *TabPage;

    void setupUi(Q3MainWindow *meshGenQt)
    {
    meshGenQt->setObjectName(QString::fromUtf8("meshGenQt"));
    widget = new QWidget(meshGenQt);
    widget->setObjectName(QString::fromUtf8("widget"));
    saveButton = new QPushButton(widget);
    saveButton->setObjectName(QString::fromUtf8("saveButton"));
    saveButton->setGeometry(QRect(11, 11, 112, 24));
    quitButton = new QPushButton(widget);
    quitButton->setObjectName(QString::fromUtf8("quitButton"));
    quitButton->setGeometry(QRect(129, 11, 112, 24));
    tabWidget = new QTabWidget(widget);
    tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
    tabWidget->setGeometry(QRect(9, 41, 580, 406));
    tab = new QWidget();
    tab->setObjectName(QString::fromUtf8("tab"));
    tabWidget->addTab(tab, QApplication::translate("meshGenQt", "General settings", 0, QApplication::UnicodeUTF8));
    tab1 = new QWidget();
    tab1->setObjectName(QString::fromUtf8("tab1"));
    tabWidget->addTab(tab1, QApplication::translate("meshGenQt", "Refinement", 0, QApplication::UnicodeUTF8));
    TabPage = new QWidget();
    TabPage->setObjectName(QString::fromUtf8("TabPage"));
    tabWidget->addTab(TabPage, QApplication::translate("meshGenQt", "Rename patches", 0, QApplication::UnicodeUTF8));
    meshGenQt->setCentralWidget(widget);

    retranslateUi(meshGenQt);

    QSize size(600, 480);
    size = size.expandedTo(meshGenQt->minimumSizeHint());
    meshGenQt->resize(size);


    QMetaObject::connectSlotsByName(meshGenQt);
    } // setupUi

    void retranslateUi(Q3MainWindow *meshGenQt)
    {
    meshGenQt->setWindowTitle(QApplication::translate("meshGenQt", "Form1", 0, QApplication::UnicodeUTF8));
    saveButton->setText(QApplication::translate("meshGenQt", "Save meshDict", 0, QApplication::UnicodeUTF8));
    quitButton->setText(QApplication::translate("meshGenQt", "Quit", 0, QApplication::UnicodeUTF8));
    tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("meshGenQt", "General settings", 0, QApplication::UnicodeUTF8));
    tabWidget->setTabText(tabWidget->indexOf(tab1), QApplication::translate("meshGenQt", "Refinement", 0, QApplication::UnicodeUTF8));
    tabWidget->setTabText(tabWidget->indexOf(TabPage), QApplication::translate("meshGenQt", "Rename patches", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(meshGenQt);
    } // retranslateUi

};

namespace Ui {
    class meshGenQt: public Ui_meshGenQt {};
} // namespace Ui

class meshGenQt : public Q3MainWindow, public Ui::meshGenQt
{
    Q_OBJECT

public:
    meshGenQt(QWidget* parent = 0, const char* name = 0, Qt::WindowFlags fl = Qt::WType_TopLevel);
    ~meshGenQt();

protected slots:
    virtual void languageChange();

};

#endif // FORM1_H
