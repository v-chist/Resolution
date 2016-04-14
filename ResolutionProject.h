#pragma once

// В файле содержится описание класса CResolutionProject - основного класса для расчета разрешения

#include "dfi.h"

#include <IC_Types.h>

#include <math.h>
#include <cmath>

#include <fstream>
#include <iostream>

#include "dfi_RasterStream.h"
#include "dfi_RasterProcs.h"
#include "dfi_raster.h"

#include "dfi_vectorObject.h"
#include "dfi_vectorObjectStream.h"

#include "dfi_style.h"
#include "VectorStyle\\IC_StyleSet.h"

#include "NoiseCalculation.h"
#include "2DVector.h"
#include "LineTarget.h"
#include "EdgeTargetNew.h"
#include "Target.h"

#include "kiss_fft.h"

class CResolutionLayer;

using namespace std;

#define R_STEP_SUCCESSFUL 1
#define STEP_UNSUCCESSFUL 0
#define L_STEP_SUCCESSFUL -1

#define STEP_X 1
#define STEP_N 2
#define STEP_R 3

#define LINES 1
#define ROWS 0


struct INTERPOLATION_POINT
{
	/*! 
	*	@struct INTERPOLATION_POINT
	*	@brief Структура для хранения точек аппроксимации резкого края.
	*/
	double x;	/*!< Положение точки относительно резкого края.*/
	double value;	/*!< Значение яркости пиксела.*/
	
	INTERPOLATION_POINT(double x_, double value_)
	{
		x		= x_;
		value	= value_;
	}
	INTERPOLATION_POINT()
	{
		x		= 0;
		value	= 0;
	}

	bool operator<(const INTERPOLATION_POINT& point)
	{
		return x < point.x;
	}
};


class CNoiseTarget
{
	/*! 
	*	@class CNoiseTarget
	*	@brief Класс для хранения результатов измерения шума.
	*/
	int m_iComboLastSelected;	

	//! Уровень сигнала.
	double m_dSignal;

	//! Уровень шума.
	double m_dDisp;				
		
public:

	void SetSignal(double dSignal);
	double GetSignal() const;
	
	void SetDisp(double dDisp);
	double GetDisp() const;

	void SetComboLastSelected(int iComboLastSelected);
	int GetComboLastSelected() const;
	
	//! Конструктор.
	CNoiseTarget();
	
	//! Конструктор.
	CNoiseTarget(double dSignal, double dDisp, int iComboLastSelected = 0);

	//! Деструктор.
	~CNoiseTarget();

	//! Оператор сравнения точек (нужен для упорядочивания точек по уровню сигнала).
	bool operator<(const CNoiseTarget& target);
};


struct CSF_PARAMETERS_SNR
{
	/*! 
	*	@struct CSF_PARAMETERS_SNR
	*	@brief Структура для хранения параметров пороговой характеристики.
	*/
	 
	double dNoise;		/*!< Отношение сигнал/шум. */
	double nTestLines;	/*!< Количество штрихов миры. */
	double dSNRLimit;	/*!< Пороговое отношение сигнал/шум. */
	
	double dMin;		/*!< Диапазон значений пикселов.*/
	double dMax;		/*!< Диапазон значений пикселов.*/ 
	
};

struct CSF_PARAMETERS 
{
	/*! 
	*	@struct CSF_PARAMETERS
	*	@brief Структура для хранения параметров пороговой характеристики.
	*/
	
	double qLimit;			/*!< Пороговое отношение сигнал/шум. */
	double m;				/*!< m = 5, 7, 9,.. (5 - для трехштриховой миры).*/
	double Dvisual;			/*!< Dза - дисперсия шума зрительного анализатора.*/
	double dGradient;		/*!< Градиент светосигнальной характеристики - параметр определяется по изображению.*/
	double dNoise;			/*!< Дисперсия шума для среднего уровня сигнала - параметр определяется по изображению.*/
	double dMaxPixValue;	/*!< Максимальное значение пиксела.*/
};

struct PROJECT_PARAMETERS
{
	/*! 
	*	@struct PROJECT_PARAMETERS
	*	@brief Структура для хранения параметров пороговой характеристики.
	*/
	int iBand;					/*!< Номер канала.*/
	int iRasterIndex;			/*!< Номер растрового слоя в документе.*/
	
	//////////////////////////////////////////////////////////////////////
	// Параметры для расчета проекции пиксела, если она не задана явно  //
	//////////////////////////////////////////////////////////////////////

	double dSatAltitude;		/*!< Высота съемки.*/
	double dFocus;				/*!< Фокусное расстояние.*/
	
	double dPixelMKMsizeX;		/*!< Размер пиксела на матрице (нм).*/
	double dPixelMKMsizeY;	
	//////////////////////////////////////////////////////////////////////

	double dPixelSizeX;			/*!< Проекция пиксела (м).*/
	double dPixelSizeY;
	
	double dContrastLimit;		/*!< Контраст тест-объекта.*/
	double dMeanPixelValue;		/*!< Среднее значение пиксела.*/
	double dLocalPixelSizeX;	/*!< Размер пиксела изображения (в локальных единицах измерения снимка).*/
	double dLocalPixelSizeY;

	////////////////////////////////////////////////////////////////////
	// Параметры выбора тест-объекта                                  //
	////////////////////////////////////////////////////////////////////
	double dDefaultTargetLength;	/*!< Размер тест-объекта по умолчанию.*/
	double dDefaultTargetWidth;	

	double dDispersePercent;		/*!< Допустимый процент дисперсии полок (процент берется от высоты резкого края).*/
	double dMinTargetWidth;			/*!< Минимальный размер тест-объекта.*/
	double dMinTargetRangePercent;	/*!< Минимальная допустимая высота тест-объекта (разница между уровнями верхней и нижней полки) в процентах от диапазона значений изображения.*/
	
	CSF_PARAMETERS_SNR CSFparamsSNR;	/*!< Cтруктура для хранения параметров пороговой характеристики.*/
};

////////////////////////////////////////////////////////////////
///// CResolutionProject
////////////////////////////////////////////////////////////////

class CResolutionProject
{
	/*! 
	*	@class CResolutionProject.
	*	@brief Класс, содержащий основную информацию о слое тест-объектов (список тест-объектов, параметры и т.д.).
	*/

	//! Указатель на объект слоя CResolutionSarLayer (слой тест-объектов).
	CResolutionLayer *m_pResolutionLayer; 

	//! Структрура, содержащая параметры слоя тест-объектов.
	PROJECT_PARAMETERS m_projectParams;

	//! Вектор тест-объектов.
	//vector <CEdgeTarget> m_listOfEdges;

	vector <CTarget*> m_listOfTargets;

	// Вектор точек шумовой характеристики.
	vector <CNoiseTarget> m_vectorNoise; 

	//! Выставление корректного знака по координате X для массива точек, аппроксимирующих резкий край.
	/*!
	* Функция выставляет для точек из resultingListOfPoints знак расстояния до резкого края таким образом,	
	*  чтобы отрицательные значения парметра "x"("расстояния") были отрицательными для нижней полки
	*	@param resultingListOfPoints - вектор точек, который будет изменен.
	*/
	void SetCorrectSign(vector<INTERPOLATION_POINT>& resultingListOfPoints) const;
	
	

	// GetFunctionalDistance возвращает значение функционала Distance, который определяется формулой:
	//             
	//            __N__                            
	//            \       /            ___      \ 2 
	// Distance =  |     |  ESF(x_i) - ESF(x_i)  |   , где 
	//            /____   \                     /
	//            i = 0
	//          ___
	//   (x_i,  ESF(x_i)) - точки, по которым производится аппроксимация (полученные по изображению), передаваемые через параметр "interpolationPoints"
	//   ESF(x)           - аналитическая аппроксимация, построенная по параметрам "params"
	//   N                - количество точек в векторе "interpolationPoints"
	//
	//сумма квадратов расстояний между измерениями ESF (по изображению) и аппроксимацией ESF с параметрами params (используется в методе наименьших квадратов)
	double GetFunctionalDistance(const vector<INTERPOLATION_POINT>& interpolationPoints, const ESF_PARAMETERS_NEW &params) const;

	double GetFunctionalDistanceLSF(const vector<INTERPOLATION_POINT>& interpolationPoints, const ESF_PARAMETERS_NEW &params, bool bBrightLine) const;

	//! Вспомогательная функция, вычисляет подыинтегральное выражение в формуле для ESF(x).
	static double Exp(const ESF_PARAMETERS_NEW& param, double z);

	//! Вспомогательная функция, вычисляет подыинтегральное выражение в формуле для LSF(x).
	static double Exp0(const ESF_PARAMETERS_NEW& param, double z);

	//! Функция, определяющая значение частотно-контрастной характеристики по известному набору значений параметров аппроксимации.
	CIC_ComplexDouble MTFValue(const ESF_PARAMETERS_NEW& param, double x) const;

	//! Вспомогательная функция для расчета значения переходной функции.
	double g(const ESF_PARAMETERS_NEW& param, double x) const;

	//! Получение минимального и максимального значения параметра x  по набору точек.
	void getMinMaxX(const vector<INTERPOLATION_POINT> &list, double& xMin, double& xMax) const;
	
	//! Функция производит один шаг итеративного процесса аппроксимации резкого края.
	/*
	*	Шаги производятся по одному из трех параметров: x0 (положение центра), R (параметр, определяющий ширину функции рассеяния линии), N (показатель степени аппроксимации).
	*	То, по какому параметру будет сделан шаг, определяется значением переменной iStepIndex. 
	* Шаг признается успешным, если значение функционала GetFunctionalDistance(..) для нового набора параметров меньше 
	* чем исходное значение функционала (dInitialDistance)
	*	@param points - набор точек, образующих переходную функцию (резкий край);
	*	@param param - параметры аппроксимации (предполагается, что параметры кроме R и x0 уже расчитаны и содержатся в структуре param).
	*	@param dInitialDistance - начальное значение функционала расстояния (входной параметр);
	*	@param dResultingDistance - результирующее значение функционала расстояния (выходной параметр);
	*	@param dStep - величина шага;
	*	@param iStepIndex - параметр определяет то, по какому параметру будет выполняться шаг:
	*
	*	• STEP_X 1 - шаг будет выполнен по параметру x0;
	*
	*	• STEP_N 2 - шаг будет выполнен по параметру N;
	*
	*	• STEP_R 3 - шаг будет выполнен по параметру R.
	*	@return Значение произведенного шага. Возможные значения: dStep, -dStep, 0. 
	* Значение 0 возвращается в том случае, если шаг по каждому из направлений по выбранной переменной не уменьшил значение функционала, определяющего точность аппроксимации.
	*/
	double DoStep(	const vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param,
								const double InitialDistance, double& resultingDistance, const double dStep, const int stepIndex) const;
	
	double DoStepLSF(	const vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, 
							const double dInitialDistance, double &dResultingDistance, const double dStep, const int iStepIndex, bool bBrightLine) const;

	//! Автоматическое определение положения резкого края.
	/*!
	* @param point - геокоординаты точки, выбранной пользователем;
	* @param R - радиус поиска (в геокоординатах);
	*/
	double ApproximateAngleDetection(const PARE &point, double Rpix) const;

	double ApproximateAngleDetectionLSF(const PARE &point, double Rpix) const;

	//! Функция определения положения резкого края.
	/*!
	* @param edgeTarget - резкий край, положение которого будет изменено;
	* @param points -  массив точек, по которым строится аппроксимация положения резкого края по методу наименьших квадратов;
	* @param nPoints - количество точек в массиве nPoints;
	* @param bDirection - параметр, характеризующий ориентацию резкого края (ROWS <=> направление резкого края ближе к вертикальному, COLUMNS <=> направление резкого края ближе к горизонтальному);
	* @param dLength - длина тест-объекта (в геокоординатах);
	* @param dWidth - ширины тест-объекта (в геокоординатах).
	*/
	void GetEdgePosition(CEdgeTargetNew *edgeTarget, const PARE *points, int nPoints, bool bDirection, double dLength, double dWidth) const;

	void GetLinePosition(CLineTarget *edgeTarget, const PARE *points, int nPoints, bool bDirection, double dLength, double dWidth) const;

	//! Функция возвращаем массив геокоординат точек, аппроксимирующих положение резкого края;
	/*!
	* @param geoRect - геокоординаты прямоугольника, внутри которого происходит поиск резкого края;
	* @param bDirection - параметр, характеризующий направление чтения изображения (LINES - считываем по строкам, ROWS - считываем по столбцам);
	* @param points -  массив точек, по которым строится аппроксимация положения резкого края по методу наименьших квадратов;
	* @param nPoints - количество точек в массиве nPoints.
	*/
	void GetLinesMaxGradPoints(const CIC_Rect3DD &geoRect, PARE *points, int nPoints, bool bDirection, bool bAbsoluteMaxCalculation = true) const;
	
	//! Функция для визуализации массива точек (заданных геокоординатами).
	void PointsVisualization(const PARE *points, int nPoints) const;

	//! Расчет начального приближения параметров аппроксимации резкого края.
	/*!
	*	@param points - набор точек, образующих переходную функцию (резкий край);
	*	@param param - параметры аппроксимации резкого края (выходной параметр).
	*/
	void GetInitialESFParameters(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param ) const;

	void InterpolatePointsLSF(vector<INTERPOLATION_POINT> &points) const;

	void InterpolatePoints(vector<INTERPOLATION_POINT> &points) const;

	void GetInitialLSFParameters(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param ) const;

	//! Расчет начальных значений параметров R и x0 аппроксимации резкого края.
	/*!
	* Предполагается, что x_min, x_max, minValue, maxValue - уже расчитаны, а param.x0 = 0.
	*	@param points - набор точек, образующих переходную функцию (резкий край);
	*	@param R - результат начальной аппроксимации (выходной параметр);
	*	@param x0 - результат начальной аппроксимации (выходной параметр);
	*	@param param - параметры аппроксимации (предполагается, что параметры кроме R и x0 уже расчитаны и содержатся в структуре param).
	*/
	void GetInitialRnX0Value(vector<INTERPOLATION_POINT> &points, double &R, double &x0, const ESF_PARAMETERS_NEW& param) const;
	
public:
	
	////////////////////////////////////////////////////////
	// Функции получения и выставления параметров класса  //
	////////////////////////////////////////////////////////
	
	PROJECT_PARAMETERS GetParams() const;			/*!< Получение структуры, содержащей параметры слоя тест-объектов.*	/

	CString GetDocName() const;						/*!< Получение имени документа.*/
	CString GetBandName() const;					/*!< Получение имени канала.*/
	CString GetRasterLayerName() const;				/*!< Получение имени растрового слоя.*/	
	CString GetResolutionLayerName() const;			/*!< Получение имени слоя тест-объекта.*/
	int GetBand() const;							/*!< Получение номера канала.*/			
	int GetRasterIndex() const;						/*!< Получение индекса растрового слоя.*/
	Cdfi_Document *GetDocument() const;				/*!< Получение указателя на документ.*/
	Cdfi_Raster *GetRaster() const;					/*!< Получение указателя на растр.*/
	CDFI *GetCdfi() const;							/*!< Получение указателя на объект класса для работы с интерфейсом.*/
	double GetContrastLimit() const;				/*!< Получение значения контраста тест-объекта.*/
	double GetPixelSizeX() const;					/*!< Получение ширины пиксела в метрах.*/
	double GetPixelSizeY() const;					/*!< Получение высоты пиксела в метрах.*/
	double GetLocalPixelSizeX() const;				/*!< Получение ширины пиксела в геокоординатах документа.*/
	double GetLocalPixelSizeY() const;				/*!< Получение высоты пиксела в геокоординатах документа.*/
	double GetDefaultTargetLength() const;			/*!< Получение длины тест-объекта по умолчанию (в пикселах).*/
	double GetDefaultTargetWidth() const;			/*!< Получение ширины тест-объекта по умолчанию (в пикселах).*/

	void SetRasterIndex(int iRasterInd);			/*!< Назначение индекса растрового слоя.*/
	void SetBand(int iBand);						/*!< Назначение номера канала.*/
	void SetResolutionLayer(CResolutionLayer *pResolutionLayer);	/*!< Назначение указателя на объект слоя CResolutionSarLayer (слой тест-объектов).*/
	CResolutionLayer *GetResolutionLayer();			/*!< Получение указателя на объект слоя CResolutionSarLayer (слой тест-объектов).*/
	void SetDefaultTargetLength(double length);		/*!< Выставление значения параметра длины резкого края по умолчанию (в пикселах).*/
	void SetDefaultTargetWidth(double width);		/*!< Выставление значения параметра ширины резкого края по умолчанию (в пикселах).*/
	void SetLocalPixelSizeX(double sizeX);			/*!< Получение ширины пиксела в геокоординатах документа.*/
	void SetLocalPixelSizeY(double sizeY);			/*!< Получение высоты пиксела в геокоординатах документа.*/

	//////////////////////////////////////////////////
	// операции по работе с вектором тест-объектов	//
	//////////////////////////////////////////////////
	CTarget *GetTarget(int i);						/*!< Получение указателя на тест-объект с индексом i.*/		
	void GetTarget(CTarget *target, int i) const;	/*!< Получение указателя на тест-объект с индексом i.*/		
	void DeleteTarget(int i);						/*!< Удаление тест-объекта с индексом i.*/
	void AddTarget(CTarget* target);				/*!< Добавление тест-объекта target в список тест-объектов.*/
	void AddEdgeTarget(const PARE &point);			/*!< Добавление тест-объекта по точке с геокоординатами point.*/
	void AddLineTarget(const PARE &point);
	void AddEdgeTarget(CEdgeTargetNew *target);
	void AddLineTarget(CLineTarget *target);

	int GetNumberOfTargets() const;					/*!< Получение количества тест-объектов.*/
	void ClearListOfTargets();						/*!< Очистить список тест-объектов.*/
	
	//////////////////////////////////////////////////////////////////
	// операции по работе с вектором точек шумовой характеристики	//
	//////////////////////////////////////////////////////////////////
	int GetNumberOfNoiseTargets();					/*!< Получение количества измерений в таблице оценки шумовой характеристики.*/
	CNoiseTarget *GetNoiseTarget(int i);			/*!< Получить результат измерения шумовой характеристики с номером i.*/
	
	//////////////////////////////////////////////////////
	// Функции для определения пригодности тест-объекта //
	//////////////////////////////////////////////////////

	//! Выставление значения флага пригодности для тест-объекта.
	bool SetValidStatus(CTarget *target) const;					
	//! Определение ширины функции рассеяния линии.
	static double ESFWidth(const ESF_PARAMETERS_NEW& param);			
	//! Расчет дисперсии полок функции рассеяния края.
	void CalculateDisp(double &dDisp, CTarget *target) const;	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//! Функция, определяющая значение функции рассеяния линии (LSF) по известному набору значений параметров аппроксимации.
	static double LSFValue(const ESF_PARAMETERS_NEW& param, double x);
	
	//! Функция, определяющая значение переходной функции (ESF) по известному набору значений параметров аппроксимации.
	static double ESFValue(const ESF_PARAMETERS_NEW& param, double x, double dDelta); 
	
	//! Функция расчета массива значений частотно-контрастной характеристики по массиву функций рассеяния линии.
	void CalculateMTF(ESF_PARAMETERS_NEW *param, int nESFplots, double dFreqMax, int N, kiss_fft_cpx *out) const;
	
	//! Функция расчета массива значений частотно-контрастной характеристики.
	void CalculateMTF(const ESF_PARAMETERS_NEW& param, double dFreqMax, int N, kiss_fft_cpx *out) const;
	
	//! Функция визуального восприятия (пороговая характеристика, функция контрастной чувствительности)
	double CSFValue(double x) const;

	//! Заполнение массива точек переходной функции (resultingListOfPoints) по тест-объекту (target).
	void FillListOfPoints(CTarget *target, vector<INTERPOLATION_POINT>& resultingListOfPoints) const;
	
	//! Уточнение значений верхней и нижней полок методом наименьших квадратов.
	/*!
	*	@param list - набор точек, образующих переходную функцию (резкий край);
	*	@param param - параметры аппроксимации резкого края (выходной параметр).
	*/
	void AdjustMinMaxValue(const vector<INTERPOLATION_POINT>& list, ESF_PARAMETERS_NEW& param) const;

	//! Расчет параметров аппроксимации резкого края.
	/*!
	*	@param target - указатель на тест объект (резкий край);
	*	@param param - параметры аппроксимации резкого края (выходной параметр).
	*/
	void GetESFapproximationParamByTarget(CEdgeTargetNew *target, ESF_PARAMETERS_NEW& param ) const;
	
	void GetLSFapproximationParamByTarget(CLineTarget *target, ESF_PARAMETERS_NEW& param ) const;

	//! Расчет параметров аппроксимации резкого края.
	/*!
	*	@param points - набор точек, образующих переходную функцию (резкий край);
	*	@param param - параметры аппроксимации резкого края (выходной параметр).
	*/
	void GetESFapproximationParamByPoints(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param/*, vector<INTERPOLATION_POINT> &outputPoints*/ ) const;

	void GetLSFapproximationParamByPoints(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, bool bBrightLine/*, vector<INTERPOLATION_POINT> &outputPoints*/, bool bFirstIteration ) const;

	//! Вычисление разрешения по указанному тест-объекту.
	void CalculateResolution(CEdgeTargetNew *target) const;

	void CalculateResolution(CLineTarget *target) const;

	//! Определение разрешения по точке (положение резкого края определяется автоматически).
	void CalculateResolutionEdge(const PARE& point) const;

	//! Автоматическое определение положения резкого края в окрестности точки, заданной параметром "point"
	void AutoEdgeDetection(const PARE &point, CEdgeTargetNew* target, int iSearchRadius) const;
	void AutoLineDetection(const PARE &point, CLineTarget* target, int iSearchRadius) const;


	//выводит результаты аппроксимации в файл
	void OutputESFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const;
	void OutputLSFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const;
	void OutputLSFPoints(const ESF_PARAMETERS_NEW& param, const vector<INTERPOLATION_POINT> &list, ofstream &resultingFile) const;
	void OutputMTFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const;
	void OutputCSFApproximationResults(ofstream &resultingFile) const;
	void OutputListOfPoints(const vector<INTERPOLATION_POINT>& list, ofstream &resultingFile) const; //выводит в файл тестовые точки
	
	//! Возвращает описывающий прямоугольник для набора тест-объектов.
	int GetBoundRect(CIC_Rect3DD& rect);

	//! Возвращает описывающий прямоугольник для выбранных тест-объектов.
	int GetBoundRectSelected(CIC_Rect3DD& rect);

	//! Функция отрисовки класса (вызывает функции отрисовки для каждого тест-объекта).
	BOOL Draw(Cdfi_View *pView, Cdfi_MemoryDC *pMemDC);

	//! Определяет, попала ли точка внутрь какого-то из тест-объектов.
	int HitTest(CPare& ptGeo, int& hitEdgeIndex);
	
	//! Выставляет флаг выделения для тест-объектов в соответствии с таблицей тест-объектов.
	void UpdateSelectionStatusFromList();

	//!	Сохраняет объект в память или читает из памяти.	
	/*!
	*	@param	pMemory - указатель на область памяти;
	*	@param	IsReading - флаг, определяющий действие с объектом ( читаем из памяти - TRUE, сохраняем в память - FALSE);
	*	@return	Размер объекта в байтах.
	*/
	long Serialize(void* pMemory, BOOL IsReading);

	//! Инициализация параметров класса.
	void Initialization(const PROJECT_PARAMETERS &params, vector <CNoiseTarget> vectorNoise, CResolutionLayer *pResolutionLayer);

	//! Конструктор.
	CResolutionProject(const PROJECT_PARAMETERS &params, CResolutionLayer *pResolutionLayer);

	//! Конструктор.
	CResolutionProject(CResolutionLayer *pResolutionLayer);

	//! Конструктор.
	CResolutionProject();

	//! Деструктор.
	~CResolutionProject();
	
	//! Оператор копирования.
	virtual void operator=(const CResolutionProject &resolutionLayer);
};



