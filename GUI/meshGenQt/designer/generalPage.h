#ifndef FORM1_H
#define FORM1_H

#include <qvariant.h>


#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QWidget>
#include <Qt3Support/Q3MimeSourceFactory>

class Ui_Form1
{
public:
    QGridLayout *gridLayout;
    QLabel *textLabel1;
    QLineEdit *lineEdit1;
    QPushButton *pushButton1;
    QLabel *textLabel2;
    QLineEdit *lineEdit2;
    QFrame *line1;
    QCheckBox *checkBox1;
    QLineEdit *lineEdit3;
    QLabel *textLabel3;
    QFrame *line2;
    QCheckBox *checkBox2;
    QLabel *textLabel4;
    QLineEdit *lineEdit4;
    QFrame *line3;
    QCheckBox *checkBox3;
    QCheckBox *checkBox4;

    void setupUi(QWidget *Form1)
    {
    Form1->setObjectName(QString::fromUtf8("Form1"));
    gridLayout = new QGridLayout(Form1);
    gridLayout->setSpacing(6);
    gridLayout->setMargin(11);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    textLabel1 = new QLabel(Form1);
    textLabel1->setObjectName(QString::fromUtf8("textLabel1"));
    textLabel1->setWordWrap(false);

    gridLayout->addWidget(textLabel1, 0, 0, 1, 1);

    lineEdit1 = new QLineEdit(Form1);
    lineEdit1->setObjectName(QString::fromUtf8("lineEdit1"));

    gridLayout->addWidget(lineEdit1, 0, 1, 1, 2);

    pushButton1 = new QPushButton(Form1);
    pushButton1->setObjectName(QString::fromUtf8("pushButton1"));

    gridLayout->addWidget(pushButton1, 0, 3, 1, 1);

    textLabel2 = new QLabel(Form1);
    textLabel2->setObjectName(QString::fromUtf8("textLabel2"));
    textLabel2->setWordWrap(false);

    gridLayout->addWidget(textLabel2, 1, 0, 1, 1);

    lineEdit2 = new QLineEdit(Form1);
    lineEdit2->setObjectName(QString::fromUtf8("lineEdit2"));

    gridLayout->addWidget(lineEdit2, 1, 1, 1, 2);

    line1 = new QFrame(Form1);
    line1->setObjectName(QString::fromUtf8("line1"));
    line1->setFrameShape(QFrame::HLine);
    line1->setFrameShadow(QFrame::Sunken);

    gridLayout->addWidget(line1, 2, 0, 1, 4);

    checkBox1 = new QCheckBox(Form1);
    checkBox1->setObjectName(QString::fromUtf8("checkBox1"));

    gridLayout->addWidget(checkBox1, 3, 0, 1, 2);

    lineEdit3 = new QLineEdit(Form1);
    lineEdit3->setObjectName(QString::fromUtf8("lineEdit3"));

    gridLayout->addWidget(lineEdit3, 4, 2, 1, 1);

    textLabel3 = new QLabel(Form1);
    textLabel3->setObjectName(QString::fromUtf8("textLabel3"));
    textLabel3->setWordWrap(false);

    gridLayout->addWidget(textLabel3, 4, 0, 1, 1);

    line2 = new QFrame(Form1);
    line2->setObjectName(QString::fromUtf8("line2"));
    line2->setFrameShape(QFrame::HLine);
    line2->setFrameShadow(QFrame::Sunken);

    gridLayout->addWidget(line2, 5, 0, 1, 4);

    checkBox2 = new QCheckBox(Form1);
    checkBox2->setObjectName(QString::fromUtf8("checkBox2"));

    gridLayout->addWidget(checkBox2, 6, 0, 1, 1);

    textLabel4 = new QLabel(Form1);
    textLabel4->setObjectName(QString::fromUtf8("textLabel4"));
    textLabel4->setWordWrap(false);

    gridLayout->addWidget(textLabel4, 7, 0, 1, 1);

    lineEdit4 = new QLineEdit(Form1);
    lineEdit4->setObjectName(QString::fromUtf8("lineEdit4"));

    gridLayout->addWidget(lineEdit4, 7, 1, 1, 2);

    line3 = new QFrame(Form1);
    line3->setObjectName(QString::fromUtf8("line3"));
    line3->setFrameShape(QFrame::HLine);
    line3->setFrameShadow(QFrame::Sunken);

    gridLayout->addWidget(line3, 8, 0, 1, 4);

    checkBox3 = new QCheckBox(Form1);
    checkBox3->setObjectName(QString::fromUtf8("checkBox3"));

    gridLayout->addWidget(checkBox3, 9, 0, 1, 3);

    checkBox4 = new QCheckBox(Form1);
    checkBox4->setObjectName(QString::fromUtf8("checkBox4"));

    gridLayout->addWidget(checkBox4, 10, 0, 1, 2);


    retranslateUi(Form1);

    QSize size(600, 480);
    size = size.expandedTo(Form1->minimumSizeHint());
    Form1->resize(size);


    QMetaObject::connectSlotsByName(Form1);
    } // setupUi

    void retranslateUi(QWidget *Form1)
    {
    Form1->setWindowTitle(QApplication::translate("Form1", "Form1", 0, QApplication::UnicodeUTF8));
    textLabel1->setText(QApplication::translate("Form1", "textLabel1", 0, QApplication::UnicodeUTF8));
    pushButton1->setText(QApplication::translate("Form1", "pushButton1", 0, QApplication::UnicodeUTF8));
    textLabel2->setText(QApplication::translate("Form1", "textLabel2", 0, QApplication::UnicodeUTF8));
    checkBox1->setText(QApplication::translate("Form1", "checkBox1", 0, QApplication::UnicodeUTF8));
    textLabel3->setText(QApplication::translate("Form1", "textLabel3", 0, QApplication::UnicodeUTF8));
    checkBox2->setText(QApplication::translate("Form1", "checkBox2", 0, QApplication::UnicodeUTF8));
    textLabel4->setText(QApplication::translate("Form1", "textLabel4", 0, QApplication::UnicodeUTF8));
    checkBox3->setText(QApplication::translate("Form1", "checkBox3", 0, QApplication::UnicodeUTF8));
    checkBox4->setText(QApplication::translate("Form1", "checkBox4", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(Form1);
    } // retranslateUi

};

namespace Ui {
    class Form1: public Ui_Form1 {};
} // namespace Ui

class Form1 : public QWidget, public Ui::Form1
{
    Q_OBJECT

public:
    Form1(QWidget* parent = 0, const char* name = 0, Qt::WindowFlags fl = 0);
    ~Form1();

protected slots:
    virtual void languageChange();

};

#endif // FORM1_H
