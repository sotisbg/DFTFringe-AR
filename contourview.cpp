/******************************************************************************
**
**  Copyright 2016 Dale Eason
**  This file is part of DFTFringe
**  is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation version 3 of the License

** DFTFringe is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with DFTFringe.  If not, see <http://www.gnu.org/licenses/>.

****************************************************************************/
#include "contourview.h"
#include "ui_contourview.h"
#include <QMenu>
#include "math.h"
#include "pixelstats.h"
#include <QSettings>
#include <QGuiApplication>
#include <QScreen>
#include <QLabel>
#include <QApplication>
contourView::contourView(QWidget *parent, ContourTools *tools) :
    QWidget(parent),
    zoomed(false), ui(new Ui::contourView), tools(tools)
{
    ui->setupUi(this);
    ui->widget->setTool(tools);
    QSettings set;
    ui->doubleSpinBox->setValue(set.value("contourRange", .100).toDouble());
    ui->fillContourCB->setChecked(set.value("contourShowFill", true).toBool());
    ui->showRuler->setChecked(set.value("contourShowRuler",false).toBool());
    setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this, &QWidget::customContextMenuRequested, this,
            &contourView::showContextMenu);
    setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    ps = new pixelStats;

    ui->LinkProfileCB->setChecked(set.value("linkProfilePlot", true).toBool());
}

contourView::~contourView()
{
    delete ui;
}
void contourView::zoom(){
    zoomed = !zoomed;
    emit zoomMe(zoomed);
}

QImage contourView::getPixstatsImage(){
    //resize(3000,2000);
        int height = QGuiApplication::primaryScreen()->geometry().height() * .75;
    QImage psImage = QImage(height, height,QImage::Format_ARGB32 );
    QPainter p3(&psImage);
    QSize originalSize = ps->size();
    ps->resize(height * .7, height);
    // The layout doesn't settle synchronously just from resize() when ps has
    // never been shown - deliver it now so imageLabel->size() below reflects
    // the size we just asked for, not a stale/default one.
    qApp->processEvents();

    // The slope error circle's label uses setScaledContents(true), which stretches
    // its (square) pixmap to fill whatever box the layout gives it - here a taller
    // box than the mirror data is wide, which is what turned the circle into an
    // ellipse in the report. Re-scale it ourselves, keeping the aspect ratio, just
    // for this capture, then put the label back the way it was.
    QLabel *imageLabel = ps->findChild<QLabel*>("image");
    QPixmap originalPixmap;
    if (imageLabel){
        originalPixmap = imageLabel->pixmap(Qt::ReturnByValue);
        if (!originalPixmap.isNull()){
            imageLabel->setScaledContents(false);
            imageLabel->setPixmap(originalPixmap.scaled(imageLabel->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
        }
    }

    ps->render(&p3);

    if (imageLabel && !originalPixmap.isNull()){
        imageLabel->setPixmap(originalPixmap);
        imageLabel->setScaledContents(true);
    }
    ps->resize(originalSize);

    return psImage;
}


void contourView::setSurface(const wavefront *wf){
    getPlot()->setSurface(wf);
    ps->setData(wf);
}

void contourView::showContextMenu(QPoint pos)
{
    // Handle global position
    QPoint globalPos = mapToGlobal(pos);

    // Create menu and insert some actions
    QMenu myMenu;
    QString txt = (zoomed)? tr("Restore to MainWindow") : tr("FullScreen");
    myMenu.addAction(txt,  this, &contourView::zoom);

    // Show context menu at handling position
    myMenu.exec(globalPos);
}

ContourPlot *contourView::getPlot(){
    return ui->widget;
}

void contourView::on_doubleSpinBox_valueChanged(double arg1)
{
    ui->widget->showContoursChanged(arg1);
    if (arg1 == 0.){
        ui->fillContourCB->setChecked(true);
        on_fillContourCB_clicked(true);
    }

}


void contourView::on_pushButton_pressed()
{
    emit showAllContours();
}

void contourView::on_histogram_clicked()
{
    ps->show();

}

void contourView::on_fillContourCB_clicked(bool checked)
{
    QSettings set;
    ui->widget->showSpectrogram(checked);

}

void contourView::updateRuler(){
    if (getPlot()->m_wf){
        QSettings settings;
        getPlot()->m_rulerPen = QPen(QColor(settings.value("ContourRulerColor", "grey").toString()));
        getPlot()->m_radialDeg = settings.value("contourRulerRadialDeg",30).toInt();
        setSurface(getPlot()->m_wf);
    }
}

void contourView::on_showRuler_clicked(bool checked)
{
    QSettings set;
    set.setValue("contourShowRuler", checked);
    if (getPlot()->m_wf)
        setSurface(getPlot()->m_wf);
}



void contourView::on_LinkProfileCB_clicked(bool checked)
{
    QSettings set;
    set.setValue("linkProfilePlot", checked);
    getPlot()->m_linkProfile = checked;
}
