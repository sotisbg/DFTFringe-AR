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
#include <qwt_plot_renderer.h>
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
    // Built from ps's own data (the histogram plot, the slope caption, the
    // native square slope-error image) instead of resizing and rendering the
    // live pixelStats widget itself. That widget is shared with the
    // interactive "Pixel Histogram" window (ps->show() elsewhere) - resizing
    // it here left it however this last left it, forcing the user to shrink
    // it back down by hand before they could use it again. Composing our own
    // image also sidesteps the label's setScaledContents(true) stretching the
    // (square) slope data into an ellipse, and any layout-timing mismatch
    // between the size we asked for and the size Qt actually laid out.
    int height = QGuiApplication::primaryScreen()->geometry().height() * .75;
    int width = height * .7;
    QImage psImage(width, height, QImage::Format_ARGB32);
    psImage.fill(Qt::white);
    QPainter p(&psImage);

    int histoHeight = height * .35;
    QwtPlotRenderer renderer;
    renderer.render(ps->histoPlot(), &p, QRect(0, 0, width, histoHeight));

    int y = histoHeight + 10;
    QFont captionFont = p.font();
    captionFont.setBold(true);
    captionFont.setPointSize(12);
    p.setFont(captionFont);
    QRect captionRect(10, y, width - 20, 60);
    p.drawText(captionRect, Qt::AlignHCenter | Qt::TextWordWrap, ps->slopeCaption());
    y += 60;

    QImage slope = ps->slopeImage();
    if (!slope.isNull()){
        int side = qMin(width - 20, height - y - 10);
        if (side > 0){
            QImage scaled = slope.scaled(side, side, Qt::KeepAspectRatio, Qt::SmoothTransformation);
            int x = (width - scaled.width()) / 2;
            p.drawImage(x, y, scaled);
            p.setPen(Qt::black);
            p.drawRect(x, y, scaled.width() - 1, scaled.height() - 1);
        }
    }

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
