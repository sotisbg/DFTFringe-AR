#include "oglrendered.h"
#include "ui_oglrendered.h"

oglRendered::oglRendered(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::oglRendered)
{
    ui->setupUi(this);

    // This widget only exists to compose the 3D surface image for the PDF
    // report, never shown interactively, so it's safe to size its captions
    // for a printed page rather than the small default UI font.
    QFont f = ui->Title->font();
    f.setPointSize(24);
    ui->Title->setFont(f);
    ui->label->setFont(f);
}

oglRendered::~oglRendered()
{
    delete ui;
}
QLabel *oglRendered::getModel(){
    return ui->model;
}
QLabel *oglRendered::getLegend(){
    return ui->legend;
}
