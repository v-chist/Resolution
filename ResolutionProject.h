#pragma once

// � ����� ���������� �������� ������ CResolutionProject - ��������� ������ ��� ������� ����������

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
	*	@brief ��������� ��� �������� ����� ������������� ������� ����.
	*/
	double x;	/*!< ��������� ����� ������������ ������� ����.*/
	double value;	/*!< �������� ������� �������.*/
	
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
	*	@brief ����� ��� �������� ����������� ��������� ����.
	*/
	int m_iComboLastSelected;	

	//! ������� �������.
	double m_dSignal;

	//! ������� ����.
	double m_dDisp;				
		
public:

	void SetSignal(double dSignal);
	double GetSignal() const;
	
	void SetDisp(double dDisp);
	double GetDisp() const;

	void SetComboLastSelected(int iComboLastSelected);
	int GetComboLastSelected() const;
	
	//! �����������.
	CNoiseTarget();
	
	//! �����������.
	CNoiseTarget(double dSignal, double dDisp, int iComboLastSelected = 0);

	//! ����������.
	~CNoiseTarget();

	//! �������� ��������� ����� (����� ��� �������������� ����� �� ������ �������).
	bool operator<(const CNoiseTarget& target);
};


struct CSF_PARAMETERS_SNR
{
	/*! 
	*	@struct CSF_PARAMETERS_SNR
	*	@brief ��������� ��� �������� ���������� ��������� ��������������.
	*/
	 
	double dNoise;		/*!< ��������� ������/���. */
	double nTestLines;	/*!< ���������� ������� ����. */
	double dSNRLimit;	/*!< ��������� ��������� ������/���. */
	
	double dMin;		/*!< �������� �������� ��������.*/
	double dMax;		/*!< �������� �������� ��������.*/ 
	
};

struct CSF_PARAMETERS 
{
	/*! 
	*	@struct CSF_PARAMETERS
	*	@brief ��������� ��� �������� ���������� ��������� ��������������.
	*/
	
	double qLimit;			/*!< ��������� ��������� ������/���. */
	double m;				/*!< m = 5, 7, 9,.. (5 - ��� ������������� ����).*/
	double Dvisual;			/*!< D�� - ��������� ���� ����������� �����������.*/
	double dGradient;		/*!< �������� ��������������� �������������� - �������� ������������ �� �����������.*/
	double dNoise;			/*!< ��������� ���� ��� �������� ������ ������� - �������� ������������ �� �����������.*/
	double dMaxPixValue;	/*!< ������������ �������� �������.*/
};

struct PROJECT_PARAMETERS
{
	/*! 
	*	@struct PROJECT_PARAMETERS
	*	@brief ��������� ��� �������� ���������� ��������� ��������������.
	*/
	int iBand;					/*!< ����� ������.*/
	int iRasterIndex;			/*!< ����� ���������� ���� � ���������.*/
	
	//////////////////////////////////////////////////////////////////////
	// ��������� ��� ������� �������� �������, ���� ��� �� ������ ����  //
	//////////////////////////////////////////////////////////////////////

	double dSatAltitude;		/*!< ������ ������.*/
	double dFocus;				/*!< �������� ����������.*/
	
	double dPixelMKMsizeX;		/*!< ������ ������� �� ������� (��).*/
	double dPixelMKMsizeY;	
	//////////////////////////////////////////////////////////////////////

	double dPixelSizeX;			/*!< �������� ������� (�).*/
	double dPixelSizeY;
	
	double dContrastLimit;		/*!< �������� ����-�������.*/
	double dMeanPixelValue;		/*!< ������� �������� �������.*/
	double dLocalPixelSizeX;	/*!< ������ ������� ����������� (� ��������� �������� ��������� ������).*/
	double dLocalPixelSizeY;

	////////////////////////////////////////////////////////////////////
	// ��������� ������ ����-�������                                  //
	////////////////////////////////////////////////////////////////////
	double dDefaultTargetLength;	/*!< ������ ����-������� �� ���������.*/
	double dDefaultTargetWidth;	

	double dDispersePercent;		/*!< ���������� ������� ��������� ����� (������� ������� �� ������ ������� ����).*/
	double dMinTargetWidth;			/*!< ����������� ������ ����-�������.*/
	double dMinTargetRangePercent;	/*!< ����������� ���������� ������ ����-������� (������� ����� �������� ������� � ������ �����) � ��������� �� ��������� �������� �����������.*/
	
	CSF_PARAMETERS_SNR CSFparamsSNR;	/*!< C�������� ��� �������� ���������� ��������� ��������������.*/
};

////////////////////////////////////////////////////////////////
///// CResolutionProject
////////////////////////////////////////////////////////////////

class CResolutionProject
{
	/*! 
	*	@class CResolutionProject.
	*	@brief �����, ���������� �������� ���������� � ���� ����-�������� (������ ����-��������, ��������� � �.�.).
	*/

	//! ��������� �� ������ ���� CResolutionSarLayer (���� ����-��������).
	CResolutionLayer *m_pResolutionLayer; 

	//! ����������, ���������� ��������� ���� ����-��������.
	PROJECT_PARAMETERS m_projectParams;

	//! ������ ����-��������.
	//vector <CEdgeTarget> m_listOfEdges;

	vector <CTarget*> m_listOfTargets;

	// ������ ����� ������� ��������������.
	vector <CNoiseTarget> m_vectorNoise; 

	//! ����������� ����������� ����� �� ���������� X ��� ������� �����, ���������������� ������ ����.
	/*!
	* ������� ���������� ��� ����� �� resultingListOfPoints ���� ���������� �� ������� ���� ����� �������,	
	*  ����� ������������� �������� �������� "x"("����������") ���� �������������� ��� ������ �����
	*	@param resultingListOfPoints - ������ �����, ������� ����� �������.
	*/
	void SetCorrectSign(vector<INTERPOLATION_POINT>& resultingListOfPoints) const;
	
	

	// GetFunctionalDistance ���������� �������� ����������� Distance, ������� ������������ ��������:
	//             
	//            __N__                            
	//            \       /            ___      \ 2 
	// Distance =  |     |  ESF(x_i) - ESF(x_i)  |   , ��� 
	//            /____   \                     /
	//            i = 0
	//          ___
	//   (x_i,  ESF(x_i)) - �����, �� ������� ������������ ������������� (���������� �� �����������), ������������ ����� �������� "interpolationPoints"
	//   ESF(x)           - ������������� �������������, ����������� �� ���������� "params"
	//   N                - ���������� ����� � ������� "interpolationPoints"
	//
	//����� ��������� ���������� ����� ����������� ESF (�� �����������) � �������������� ESF � ����������� params (������������ � ������ ���������� ���������)
	double GetFunctionalDistance(const vector<INTERPOLATION_POINT>& interpolationPoints, const ESF_PARAMETERS_NEW &params) const;

	double GetFunctionalDistanceLSF(const vector<INTERPOLATION_POINT>& interpolationPoints, const ESF_PARAMETERS_NEW &params, bool bBrightLine) const;

	//! ��������������� �������, ��������� ���������������� ��������� � ������� ��� ESF(x).
	static double Exp(const ESF_PARAMETERS_NEW& param, double z);

	//! ��������������� �������, ��������� ���������������� ��������� � ������� ��� LSF(x).
	static double Exp0(const ESF_PARAMETERS_NEW& param, double z);

	//! �������, ������������ �������� ��������-����������� �������������� �� ���������� ������ �������� ���������� �������������.
	CIC_ComplexDouble MTFValue(const ESF_PARAMETERS_NEW& param, double x) const;

	//! ��������������� ������� ��� ������� �������� ���������� �������.
	double g(const ESF_PARAMETERS_NEW& param, double x) const;

	//! ��������� ������������ � ������������� �������� ��������� x  �� ������ �����.
	void getMinMaxX(const vector<INTERPOLATION_POINT> &list, double& xMin, double& xMax) const;
	
	//! ������� ���������� ���� ��� ������������ �������� ������������� ������� ����.
	/*
	*	���� ������������ �� ������ �� ���� ����������: x0 (��������� ������), R (��������, ������������ ������ ������� ��������� �����), N (���������� ������� �������������).
	*	��, �� ������ ��������� ����� ������ ���, ������������ ��������� ���������� iStepIndex. 
	* ��� ���������� ��������, ���� �������� ����������� GetFunctionalDistance(..) ��� ������ ������ ���������� ������ 
	* ��� �������� �������� ����������� (dInitialDistance)
	*	@param points - ����� �����, ���������� ���������� ������� (������ ����);
	*	@param param - ��������� ������������� (��������������, ��� ��������� ����� R � x0 ��� ��������� � ���������� � ��������� param).
	*	@param dInitialDistance - ��������� �������� ����������� ���������� (������� ��������);
	*	@param dResultingDistance - �������������� �������� ����������� ���������� (�������� ��������);
	*	@param dStep - �������� ����;
	*	@param iStepIndex - �������� ���������� ��, �� ������ ��������� ����� ����������� ���:
	*
	*	� STEP_X 1 - ��� ����� �������� �� ��������� x0;
	*
	*	� STEP_N 2 - ��� ����� �������� �� ��������� N;
	*
	*	� STEP_R 3 - ��� ����� �������� �� ��������� R.
	*	@return �������� �������������� ����. ��������� ��������: dStep, -dStep, 0. 
	* �������� 0 ������������ � ��� ������, ���� ��� �� ������� �� ����������� �� ��������� ���������� �� �������� �������� �����������, ������������� �������� �������������.
	*/
	double DoStep(	const vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param,
								const double InitialDistance, double& resultingDistance, const double dStep, const int stepIndex) const;
	
	double DoStepLSF(	const vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, 
							const double dInitialDistance, double &dResultingDistance, const double dStep, const int iStepIndex, bool bBrightLine) const;

	//! �������������� ����������� ��������� ������� ����.
	/*!
	* @param point - ������������� �����, ��������� �������������;
	* @param R - ������ ������ (� ��������������);
	*/
	double ApproximateAngleDetection(const PARE &point, double Rpix) const;

	double ApproximateAngleDetectionLSF(const PARE &point, double Rpix) const;

	//! ������� ����������� ��������� ������� ����.
	/*!
	* @param edgeTarget - ������ ����, ��������� �������� ����� ��������;
	* @param points -  ������ �����, �� ������� �������� ������������� ��������� ������� ���� �� ������ ���������� ���������;
	* @param nPoints - ���������� ����� � ������� nPoints;
	* @param bDirection - ��������, ��������������� ���������� ������� ���� (ROWS <=> ����������� ������� ���� ����� � �������������, COLUMNS <=> ����������� ������� ���� ����� � ���������������);
	* @param dLength - ����� ����-������� (� ��������������);
	* @param dWidth - ������ ����-������� (� ��������������).
	*/
	void GetEdgePosition(CEdgeTargetNew *edgeTarget, const PARE *points, int nPoints, bool bDirection, double dLength, double dWidth) const;

	void GetLinePosition(CLineTarget *edgeTarget, const PARE *points, int nPoints, bool bDirection, double dLength, double dWidth) const;

	//! ������� ���������� ������ ������������ �����, ���������������� ��������� ������� ����;
	/*!
	* @param geoRect - ������������� ��������������, ������ �������� ���������� ����� ������� ����;
	* @param bDirection - ��������, ��������������� ����������� ������ ����������� (LINES - ��������� �� �������, ROWS - ��������� �� ��������);
	* @param points -  ������ �����, �� ������� �������� ������������� ��������� ������� ���� �� ������ ���������� ���������;
	* @param nPoints - ���������� ����� � ������� nPoints.
	*/
	void GetLinesMaxGradPoints(const CIC_Rect3DD &geoRect, PARE *points, int nPoints, bool bDirection, bool bAbsoluteMaxCalculation = true) const;
	
	//! ������� ��� ������������ ������� ����� (�������� ���������������).
	void PointsVisualization(const PARE *points, int nPoints) const;

	//! ������ ���������� ����������� ���������� ������������� ������� ����.
	/*!
	*	@param points - ����� �����, ���������� ���������� ������� (������ ����);
	*	@param param - ��������� ������������� ������� ���� (�������� ��������).
	*/
	void GetInitialESFParameters(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param ) const;

	void InterpolatePointsLSF(vector<INTERPOLATION_POINT> &points) const;

	void InterpolatePoints(vector<INTERPOLATION_POINT> &points) const;

	void GetInitialLSFParameters(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param ) const;

	//! ������ ��������� �������� ���������� R � x0 ������������� ������� ����.
	/*!
	* ��������������, ��� x_min, x_max, minValue, maxValue - ��� ���������, � param.x0 = 0.
	*	@param points - ����� �����, ���������� ���������� ������� (������ ����);
	*	@param R - ��������� ��������� ������������� (�������� ��������);
	*	@param x0 - ��������� ��������� ������������� (�������� ��������);
	*	@param param - ��������� ������������� (��������������, ��� ��������� ����� R � x0 ��� ��������� � ���������� � ��������� param).
	*/
	void GetInitialRnX0Value(vector<INTERPOLATION_POINT> &points, double &R, double &x0, const ESF_PARAMETERS_NEW& param) const;
	
public:
	
	////////////////////////////////////////////////////////
	// ������� ��������� � ����������� ���������� ������  //
	////////////////////////////////////////////////////////
	
	PROJECT_PARAMETERS GetParams() const;			/*!< ��������� ���������, ���������� ��������� ���� ����-��������.*	/

	CString GetDocName() const;						/*!< ��������� ����� ���������.*/
	CString GetBandName() const;					/*!< ��������� ����� ������.*/
	CString GetRasterLayerName() const;				/*!< ��������� ����� ���������� ����.*/	
	CString GetResolutionLayerName() const;			/*!< ��������� ����� ���� ����-�������.*/
	int GetBand() const;							/*!< ��������� ������ ������.*/			
	int GetRasterIndex() const;						/*!< ��������� ������� ���������� ����.*/
	Cdfi_Document *GetDocument() const;				/*!< ��������� ��������� �� ��������.*/
	Cdfi_Raster *GetRaster() const;					/*!< ��������� ��������� �� �����.*/
	CDFI *GetCdfi() const;							/*!< ��������� ��������� �� ������ ������ ��� ������ � �����������.*/
	double GetContrastLimit() const;				/*!< ��������� �������� ��������� ����-�������.*/
	double GetPixelSizeX() const;					/*!< ��������� ������ ������� � ������.*/
	double GetPixelSizeY() const;					/*!< ��������� ������ ������� � ������.*/
	double GetLocalPixelSizeX() const;				/*!< ��������� ������ ������� � �������������� ���������.*/
	double GetLocalPixelSizeY() const;				/*!< ��������� ������ ������� � �������������� ���������.*/
	double GetDefaultTargetLength() const;			/*!< ��������� ����� ����-������� �� ��������� (� ��������).*/
	double GetDefaultTargetWidth() const;			/*!< ��������� ������ ����-������� �� ��������� (� ��������).*/

	void SetRasterIndex(int iRasterInd);			/*!< ���������� ������� ���������� ����.*/
	void SetBand(int iBand);						/*!< ���������� ������ ������.*/
	void SetResolutionLayer(CResolutionLayer *pResolutionLayer);	/*!< ���������� ��������� �� ������ ���� CResolutionSarLayer (���� ����-��������).*/
	CResolutionLayer *GetResolutionLayer();			/*!< ��������� ��������� �� ������ ���� CResolutionSarLayer (���� ����-��������).*/
	void SetDefaultTargetLength(double length);		/*!< ����������� �������� ��������� ����� ������� ���� �� ��������� (� ��������).*/
	void SetDefaultTargetWidth(double width);		/*!< ����������� �������� ��������� ������ ������� ���� �� ��������� (� ��������).*/
	void SetLocalPixelSizeX(double sizeX);			/*!< ��������� ������ ������� � �������������� ���������.*/
	void SetLocalPixelSizeY(double sizeY);			/*!< ��������� ������ ������� � �������������� ���������.*/

	//////////////////////////////////////////////////
	// �������� �� ������ � �������� ����-��������	//
	//////////////////////////////////////////////////
	CTarget *GetTarget(int i);						/*!< ��������� ��������� �� ����-������ � �������� i.*/		
	void GetTarget(CTarget *target, int i) const;	/*!< ��������� ��������� �� ����-������ � �������� i.*/		
	void DeleteTarget(int i);						/*!< �������� ����-������� � �������� i.*/
	void AddTarget(CTarget* target);				/*!< ���������� ����-������� target � ������ ����-��������.*/
	void AddEdgeTarget(const PARE &point);			/*!< ���������� ����-������� �� ����� � ��������������� point.*/
	void AddLineTarget(const PARE &point);
	void AddEdgeTarget(CEdgeTargetNew *target);
	void AddLineTarget(CLineTarget *target);

	int GetNumberOfTargets() const;					/*!< ��������� ���������� ����-��������.*/
	void ClearListOfTargets();						/*!< �������� ������ ����-��������.*/
	
	//////////////////////////////////////////////////////////////////
	// �������� �� ������ � �������� ����� ������� ��������������	//
	//////////////////////////////////////////////////////////////////
	int GetNumberOfNoiseTargets();					/*!< ��������� ���������� ��������� � ������� ������ ������� ��������������.*/
	CNoiseTarget *GetNoiseTarget(int i);			/*!< �������� ��������� ��������� ������� �������������� � ������� i.*/
	
	//////////////////////////////////////////////////////
	// ������� ��� ����������� ����������� ����-������� //
	//////////////////////////////////////////////////////

	//! ����������� �������� ����� ����������� ��� ����-�������.
	bool SetValidStatus(CTarget *target) const;					
	//! ����������� ������ ������� ��������� �����.
	static double ESFWidth(const ESF_PARAMETERS_NEW& param);			
	//! ������ ��������� ����� ������� ��������� ����.
	void CalculateDisp(double &dDisp, CTarget *target) const;	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//! �������, ������������ �������� ������� ��������� ����� (LSF) �� ���������� ������ �������� ���������� �������������.
	static double LSFValue(const ESF_PARAMETERS_NEW& param, double x);
	
	//! �������, ������������ �������� ���������� ������� (ESF) �� ���������� ������ �������� ���������� �������������.
	static double ESFValue(const ESF_PARAMETERS_NEW& param, double x, double dDelta); 
	
	//! ������� ������� ������� �������� ��������-����������� �������������� �� ������� ������� ��������� �����.
	void CalculateMTF(ESF_PARAMETERS_NEW *param, int nESFplots, double dFreqMax, int N, kiss_fft_cpx *out) const;
	
	//! ������� ������� ������� �������� ��������-����������� ��������������.
	void CalculateMTF(const ESF_PARAMETERS_NEW& param, double dFreqMax, int N, kiss_fft_cpx *out) const;
	
	//! ������� ����������� ���������� (��������� ��������������, ������� ����������� ����������������)
	double CSFValue(double x) const;

	//! ���������� ������� ����� ���������� ������� (resultingListOfPoints) �� ����-������� (target).
	void FillListOfPoints(CTarget *target, vector<INTERPOLATION_POINT>& resultingListOfPoints) const;
	
	//! ��������� �������� ������� � ������ ����� ������� ���������� ���������.
	/*!
	*	@param list - ����� �����, ���������� ���������� ������� (������ ����);
	*	@param param - ��������� ������������� ������� ���� (�������� ��������).
	*/
	void AdjustMinMaxValue(const vector<INTERPOLATION_POINT>& list, ESF_PARAMETERS_NEW& param) const;

	//! ������ ���������� ������������� ������� ����.
	/*!
	*	@param target - ��������� �� ���� ������ (������ ����);
	*	@param param - ��������� ������������� ������� ���� (�������� ��������).
	*/
	void GetESFapproximationParamByTarget(CEdgeTargetNew *target, ESF_PARAMETERS_NEW& param ) const;
	
	void GetLSFapproximationParamByTarget(CLineTarget *target, ESF_PARAMETERS_NEW& param ) const;

	//! ������ ���������� ������������� ������� ����.
	/*!
	*	@param points - ����� �����, ���������� ���������� ������� (������ ����);
	*	@param param - ��������� ������������� ������� ���� (�������� ��������).
	*/
	void GetESFapproximationParamByPoints(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param/*, vector<INTERPOLATION_POINT> &outputPoints*/ ) const;

	void GetLSFapproximationParamByPoints(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, bool bBrightLine/*, vector<INTERPOLATION_POINT> &outputPoints*/, bool bFirstIteration ) const;

	//! ���������� ���������� �� ���������� ����-�������.
	void CalculateResolution(CEdgeTargetNew *target) const;

	void CalculateResolution(CLineTarget *target) const;

	//! ����������� ���������� �� ����� (��������� ������� ���� ������������ �������������).
	void CalculateResolutionEdge(const PARE& point) const;

	//! �������������� ����������� ��������� ������� ���� � ����������� �����, �������� ���������� "point"
	void AutoEdgeDetection(const PARE &point, CEdgeTargetNew* target, int iSearchRadius) const;
	void AutoLineDetection(const PARE &point, CLineTarget* target, int iSearchRadius) const;


	//������� ���������� ������������� � ����
	void OutputESFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const;
	void OutputLSFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const;
	void OutputLSFPoints(const ESF_PARAMETERS_NEW& param, const vector<INTERPOLATION_POINT> &list, ofstream &resultingFile) const;
	void OutputMTFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const;
	void OutputCSFApproximationResults(ofstream &resultingFile) const;
	void OutputListOfPoints(const vector<INTERPOLATION_POINT>& list, ofstream &resultingFile) const; //������� � ���� �������� �����
	
	//! ���������� ����������� ������������� ��� ������ ����-��������.
	int GetBoundRect(CIC_Rect3DD& rect);

	//! ���������� ����������� ������������� ��� ��������� ����-��������.
	int GetBoundRectSelected(CIC_Rect3DD& rect);

	//! ������� ��������� ������ (�������� ������� ��������� ��� ������� ����-�������).
	BOOL Draw(Cdfi_View *pView, Cdfi_MemoryDC *pMemDC);

	//! ����������, ������ �� ����� ������ ������-�� �� ����-��������.
	int HitTest(CPare& ptGeo, int& hitEdgeIndex);
	
	//! ���������� ���� ��������� ��� ����-�������� � ������������ � �������� ����-��������.
	void UpdateSelectionStatusFromList();

	//!	��������� ������ � ������ ��� ������ �� ������.	
	/*!
	*	@param	pMemory - ��������� �� ������� ������;
	*	@param	IsReading - ����, ������������ �������� � �������� ( ������ �� ������ - TRUE, ��������� � ������ - FALSE);
	*	@return	������ ������� � ������.
	*/
	long Serialize(void* pMemory, BOOL IsReading);

	//! ������������� ���������� ������.
	void Initialization(const PROJECT_PARAMETERS &params, vector <CNoiseTarget> vectorNoise, CResolutionLayer *pResolutionLayer);

	//! �����������.
	CResolutionProject(const PROJECT_PARAMETERS &params, CResolutionLayer *pResolutionLayer);

	//! �����������.
	CResolutionProject(CResolutionLayer *pResolutionLayer);

	//! �����������.
	CResolutionProject();

	//! ����������.
	~CResolutionProject();
	
	//! �������� �����������.
	virtual void operator=(const CResolutionProject &resolutionLayer);
};



