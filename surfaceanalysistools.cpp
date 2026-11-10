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
#include "surfaceanalysistools.h"
#include "ui_surfaceanalysistools.h"
#include <QAction>
#include <QDebug>
#include <QMenu>
#include <QSettings>
#include <QShortcut>
#include <QPainter>
#include <QPixmap>
#include <QtAlgorithms>
#include <QLineEdit>
#include <QSignalBlocker>
#include "mirrordlg.h"
#include <QMessageBox>
surfaceAnalysisTools *surfaceAnalysisTools::m_Instance = NULL;

surfaceAnalysisTools * surfaceAnalysisTools::get_Instance(QWidget *parent){
    if  (m_Instance == NULL){
        m_Instance = new surfaceAnalysisTools(parent);
    }
    return m_Instance;
}

surfaceAnalysisTools::surfaceAnalysisTools(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::surfaceAnalysisTools),lastCurrentItem(-1)
{
    ui->setupUi(this);
    ui->wavefrontList->setIconSize(QSize(34, 16));   // check mark + status dot, see refreshIcon()
    //ui->deleteWave->setIcon(style()->standardIcon(QStyle::SP_TrashIcon));
    m_edgeMaskLabel = ui->outsideEdgeMaskLabel;
    m_centerMaskLabel = ui->insideEdgeMaskLabel;
    m_blurrRadius = ui->blurMm;
    QSettings settings;
    bool ch = settings.value("GBlur", true).toBool();
    double val = settings.value("GBValue", 20).toInt();
    ui->surfaceSmoothGausianBlurr->setValue(val);
    ui->blurCB->setCheckState((ch) ? Qt::Checked : Qt::Unchecked);
    m_useDefocus = false;
    m_defocus = 0.;
    connect(&m_defocusTimer, &QTimer::timeout, this, &surfaceAnalysisTools::defocusTimerDone);
    connect(ui->wavefrontList->itemDelegate(), &QAbstractItemDelegate::closeEditor, this,
            &surfaceAnalysisTools::ListWidgetEditEnd);
    ui->wavefrontList->setContextMenuPolicy(Qt::CustomContextMenu);
    QShortcut* delShortcut = new QShortcut(QKeySequence::Delete, this);
    connect(delShortcut, &QShortcut::activated, this, &surfaceAnalysisTools::on_deleteWave_clicked);
}

void surfaceAnalysisTools::enableControls(bool flag){
    setEnabled(flag);
}


void surfaceAnalysisTools::setBlurText(const QString &txt){
    ui->blurMm->setText(txt);
}

void surfaceAnalysisTools::setGaussianControls(bool enabled, double value){
    const QSignalBlocker blurToggleBlocker(ui->blurCB);
    const QSignalBlocker blurValueBlocker(ui->surfaceSmoothGausianBlurr);
    ui->blurCB->setChecked(enabled);
    ui->surfaceSmoothGausianBlurr->setValue(value);
    ui->surfaceSmoothGausianBlurr->setEnabled(enabled);
}

void surfaceAnalysisTools::setGaussianStateText(const QString & str){
    ui->lblGaussState->setText(str);
}


void surfaceAnalysisTools::addWaveFront(const QString &name){

    QStringList list = name.split('/');
    QString shorter = name;
    if (list.size() > 1)
        shorter = list[list.size()-2] + "/" + list[list.size()-1];//fontMetrics().elidedText(name,Qt::ElideLeft,w);

    ui->wavefrontList->addItem(shorter);
    currentNdxChanged(ui->wavefrontList->count()-1);
    ui->wavefrontList->item(lastCurrentItem)->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
}

void surfaceAnalysisTools::removeWaveFront(const QString &){

}

surfaceAnalysisTools::~surfaceAnalysisTools()
{
    delete ui;
}

void surfaceAnalysisTools::on_blurCB_clicked(bool checked)
{
    ui->surfaceSmoothGausianBlurr->setEnabled(checked);
    emit surfaceSmoothGBEnabled(checked);
}

void surfaceAnalysisTools::on_spinBox_valueChanged(int arg1)
{
    emit outsideMaskValue(arg1);
}

void surfaceAnalysisTools::on_spinBox_2_valueChanged(int arg1)
{
    emit centerMaskValue(arg1);
}

// Icon of a row = [check mark if it is the current wave front][yellow/red dot if
// Accuracy found it deviating/far out]; both slots are always 16 px wide so the
// names stay aligned.  The flag lives in the item (Qt::UserRole + 1).
void surfaceAnalysisTools::refreshIcon(int row){
    QListWidgetItem *item = ui->wavefrontList->item(row);
    if (!item)
        return;
    int flag = item->data(Qt::UserRole + 1).toInt();
    QPixmap pm(34, 16);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setRenderHint(QPainter::Antialiasing);
    if (row == lastCurrentItem)
        p.drawPixmap(0, 0, style()->standardIcon(QStyle::SP_DialogYesButton).pixmap(16, 16));
    if (flag == 1 || flag == 2){
        p.setPen(QPen(flag == 2 ? QColor(150, 0, 0) : QColor(160, 110, 0), 1));
        p.setBrush(flag == 2 ? QColor(220, 30, 30) : QColor(250, 200, 0));
        p.drawEllipse(QRectF(19, 3, 10, 10));
    }
    p.end();
    item->setIcon(QIcon(pm));
}

// Marks the wave fronts Accuracy found deviating (yellow dot) or far out (red
// dot) and selects the far-out ones, so Delete removes them.
void surfaceAnalysisTools::markWaveFronts(const QList<int> &flags){
    ui->wavefrontList->clearSelection();
    for (int i = 0; i < ui->wavefrontList->count(); ++i){
        int f = i < flags.size() ? flags[i] : 0;
        ui->wavefrontList->item(i)->setData(Qt::UserRole + 1, f);
        refreshIcon(i);
        if (f == 2)
            ui->wavefrontList->item(i)->setSelected(true);
    }
    ui->wavefrontList->setFocus();
}

// this is responsible of giving user feedback about which wft is currently displayed
void surfaceAnalysisTools::currentNdxChanged(int ndx){
    if (ui->wavefrontList->count() == 0)
        return;

    int prev = lastCurrentItem;
    lastCurrentItem = -1;
    refreshIcon(prev);   // the old current one loses the check mark
    ui->wavefrontList->setCurrentRow(ndx,QItemSelectionModel::Current);

    lastCurrentItem = ui->wavefrontList->currentRow();
    refreshIcon(lastCurrentItem);
}

void surfaceAnalysisTools::deleteWaveFront(int i){
    if(lastCurrentItem > i) {
        lastCurrentItem--;
    }
    QListWidgetItem* item = ui->wavefrontList->takeItem(i);
    delete item;
}

// this is called when pressing enter on a wavefront selected with arrow key
void surfaceAnalysisTools::on_wavefrontList_activated(const QModelIndex &index)
{
    on_wavefrontList_clicked(index);
}
void surfaceAnalysisTools::on_wavefrontList_clicked(const QModelIndex &index)
{
    /*
    QModelIndexList indexes = ui->wavefrontList->selectionModel()->selectedIndexes();
    foreach(QModelIndex index, indexes)
    {
        emit waveFrontClicked(index.row());
    }
            */
    // only inform about the clicked item in the list for efficency.
    //That way surface manager does not spend time updating things not shown.
    emit waveFrontClicked(index.row());
}

// this is "right click" on windows
void surfaceAnalysisTools::on_wavefrontList_customContextMenuRequested(QPoint pos)
{
    QListWidgetItem *item = ui->wavefrontList->itemAt(pos);
    if (item == 0)
        return;

    QMenu contextMenu(this);

    // Create "Delete" action
    QAction *deleteAction = new QAction("Delete", this);
    connect(deleteAction, &QAction::triggered, this, [this]() {
        on_deleteWave_clicked();
    });
    deleteAction->setShortcut(QKeySequence::Delete);
    contextMenu.addAction(deleteAction);

    // Create "Invert" action
    QAction *invertAction = new QAction("Invert", this);
    connect(invertAction, &QAction::triggered, this, [this]() {
        on_InvertPB_pressed();
    });
    contextMenu.addAction(invertAction);

    // Create "Rename" action
    QAction *renameAction = new QAction("Rename", this);
    connect(renameAction, &QAction::triggered, this, [this, item]() {
        ui->wavefrontList->editItem(item);
    });
    renameAction->setShortcut(QKeySequence(Qt::Key_F2));
    contextMenu.addAction(renameAction);

    // Create "Rotate" action
    QAction *rotateAction = new QAction("Rotate", this);
    connect(rotateAction, &QAction::triggered, this, [this]() {
        on_transformsPB_clicked();
    });
    contextMenu.addAction(rotateAction);

    // Show the context menu at the cursor position
    contextMenu.exec(ui->wavefrontList->mapToGlobal(pos));
}

QList<int> surfaceAnalysisTools::SelectedWaveFronts(){
    QModelIndexList indexes = ui->wavefrontList->selectionModel()->selectedIndexes();

    QList<int> indexList;
    foreach(QModelIndex index, indexes)
    {
        indexList << index.row();
    }
    return indexList;
}

void surfaceAnalysisTools::on_deleteWave_clicked()
{
    QList<int> indexList = SelectedWaveFronts();
    std::sort(indexList.begin(), indexList.end(), std::greater<int>());
    emit deleteTheseWaveFronts(indexList);
}




void surfaceAnalysisTools::on_transformsPB_clicked()
{
    emit doxform(SelectedWaveFronts());
}

void surfaceAnalysisTools::on_averagePB_clicked()
{
    emit average(SelectedWaveFronts());
}

void surfaceAnalysisTools::on_accuracyPB_clicked()
{
    emit accuracy(SelectedWaveFronts());
}
void surfaceAnalysisTools::select(int item){
    ui->wavefrontList->selectionModel()->reset();
    if (ui->wavefrontList->count() > 0)
        ui->wavefrontList->item(item)->setSelected(true);
}

void surfaceAnalysisTools::on_SelectAllPB_clicked()
{
    ui->wavefrontList->selectAll();
}

void surfaceAnalysisTools::on_SelectNonePB_clicked()
{
    ui->wavefrontList->selectionModel()->reset();
}
void surfaceAnalysisTools::nameChangedN(int ndx, const QString &newname){
    ui->wavefrontList->item(ndx)->setText(newname);
}

void surfaceAnalysisTools::nameChanged(const QString &old, const QString &newname){

    // first look for full length name match
    QList<QListWidgetItem *> ql = ui->wavefrontList->findItems(old,Qt::MatchExactly);
    if (ql.size() == 1){
        ql[0]->setText(newname);
        return;
    }

    QStringList list = old.split('/');
    QString shorter = old;

    if (list.size() > 1)
        shorter = list[list.size()-2] + "/" + list[list.size()-1];

    ql = ui->wavefrontList->findItems(shorter,Qt::MatchEndsWith);
    if (ql.size()) {

        list = newname.split('/');
        QString shorter = newname;

        if (list.size() > 1)
            shorter = list[list.size()-2] + "/" + list[list.size()-1];
        ql[0]->setText(shorter);
    }
}



void surfaceAnalysisTools::defocusTimerDone(){
    m_defocusTimer.stop();
    emit defocusChanged();
}

#include "defocusdlg.h"
void surfaceAnalysisTools::defocusControlChanged(double val){
    m_defocus = val;
    emit defocusChanged();

}


void surfaceAnalysisTools::on_InvertPB_pressed()
{
    emit invert(SelectedWaveFronts());
}

// this is triggered on renaming/editing end
void surfaceAnalysisTools::ListWidgetEditEnd(QWidget* editor, QAbstractItemDelegate::EndEditHint ){
    const QString NewValue = reinterpret_cast<QLineEdit*>(editor)->text();
    QModelIndexList indexList = ui->wavefrontList->selectionModel()->selectedIndexes();
    // because we don't know exactly which wft has been renamed in case of multiple selection, we must update name for all
    for (const QModelIndex &index : indexList) {
        const int row = index.row();
        // update name only for thoseactually matching new name (that have actually be renamed)
        if(NewValue == ui->wavefrontList->item(row)->text()){
            emit wftNameChanged(index.row(), ui->wavefrontList->item(row)->text());
            //force complete redisplay to reload if rename was on currently displayed wft
            if(index.row() == lastCurrentItem) {
                on_wavefrontList_clicked(index);
            }
            // we cannot break; the loop because there might be duplicate names (rare but can happen)
        }
    }
}

void surfaceAnalysisTools::on_pushButton_clicked()
{
    emit updateSelected();
}


void surfaceAnalysisTools::on_filterPB_clicked()
{
    emit filterWavefronts();
}



void surfaceAnalysisTools::on_surfaceSmoothGausianBlurr_valueChanged(double arg1)
{
   emit     surfaceSmoothGBValue(arg1);
}






void surfaceAnalysisTools::on_defocus_clicked()
{
    if (ui->wavefrontList->count() == 0){
        return;
    }
    m_useDefocus = true;
    emit defocusSetup();

    defocusDlg *dlg = new defocusDlg();
    connect(dlg, &defocusDlg::defocus, this, &surfaceAnalysisTools::defocusControlChanged);
    connect(dlg, &QDialog::finished, this, &surfaceAnalysisTools::closeDefocus );
    dlg->show();

return;
}

void surfaceAnalysisTools::closeDefocus(int /*result*/){
    m_useDefocus = false;
    m_defocus = 0.;
    emit defocusChanged();
}


void surfaceAnalysisTools::deleteLast(){
    if (ui->wavefrontList->count() > 0) {
        // Calculate the index of the last item
        int lastRow = ui->wavefrontList->count() - 1;

        // Remove the item from the widget (returns a pointer to the item)
        QListWidgetItem *item = ui->wavefrontList->takeItem(lastRow);

        // Explicitly delete the item to prevent memory leaks
        delete item;
    }
}
