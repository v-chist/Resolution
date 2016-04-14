
#include "stdafx.h"

#include "Resolution.h"

#include "ResolutionProject.h"
#include "ResolutionLayer.h"
//#include "NewPrjPageParameters_SNR.h"

#include "dfi_layer.h"
#include "fft.h"
#include <IC_Types.h>
#include "kiss_fft.h"

////////////////////////////////////////////////////////////
// CResolutionProject - основной класс для расчета разрешения
////////////////////////////////////////////////////////////

CResolutionProject::CResolutionProject(const PROJECT_PARAMETERS &params, CResolutionLayer *pResolutionLayer)
{
	m_pResolutionLayer = pResolutionLayer;
	m_projectParams = params;
}

CResolutionProject::CResolutionProject(CResolutionLayer *pResolutionLayer)
{
	m_pResolutionLayer = pResolutionLayer;
}

CResolutionProject::CResolutionProject()
{
	m_pResolutionLayer = 0;
}

void CResolutionProject::operator=(const CResolutionProject &resolutionProject)
{
	m_pResolutionLayer = resolutionProject.m_pResolutionLayer;

	m_projectParams = resolutionProject.m_projectParams;

	m_projectParams.dPixelSizeX = resolutionProject.m_projectParams.dPixelSizeX;
	m_projectParams.dPixelSizeY = resolutionProject.m_projectParams.dPixelSizeY;

	m_listOfTargets.clear();

	for(int i = 0; i < resolutionProject.GetNumberOfTargets(); i++)
		m_listOfTargets.push_back(resolutionProject.m_listOfTargets[i]->MakeTargetCopy());

	//Копируем список измерений шума
	m_vectorNoise.clear();
	for(size_t i = 0; i < resolutionProject.m_vectorNoise.size(); i++)
		m_vectorNoise.push_back(resolutionProject.m_vectorNoise[i]);
}

void CResolutionProject::Initialization(const PROJECT_PARAMETERS &params, vector <CNoiseTarget> vectorNoise, CResolutionLayer *pResolutionLayer)
{
	m_projectParams = params;
	m_pResolutionLayer = pResolutionLayer;
	m_vectorNoise.clear();
	int size = vectorNoise.size();
	for(int i = 0; i < size; i++)
		m_vectorNoise.push_back(vectorNoise[i]);
}

CResolutionProject::~CResolutionProject()
{
	ClearListOfTargets();
}

CTarget* CResolutionProject::GetTarget(int i) 
{
	int iSize = m_listOfTargets.size();
	if(i >= iSize)
		return NULL;

	return m_listOfTargets[i]; 
}

void CResolutionProject::GetTarget(CTarget *target, int i) const
{
	if(target->GetType() == LINE_TARGET)
	{
		CLineTarget *lineTarget = (CLineTarget *)m_listOfTargets[i];
		
		PARE point1(lineTarget->GetPoint1());
		PARE point2(lineTarget->GetPoint2());
		((CLineTarget *)target)->Initialize( point1, point2, lineTarget->GetWidth1(), lineTarget->GetWidth2(), m_projectParams.dLocalPixelSizeX, m_projectParams.dLocalPixelSizeY);
	}
	else if(target->GetType() == EDGE_TARGET)
	{
		CEdgeTargetNew *edgeTarget = (CEdgeTargetNew *)m_listOfTargets[i];
		
		PARE point1(edgeTarget->GetPoint1());
		PARE point2(edgeTarget->GetPoint2());
		((CEdgeTargetNew *)target)->Initialize( point1, point2, edgeTarget->GetWidth1(), edgeTarget->GetWidth2(), m_projectParams.dLocalPixelSizeX, m_projectParams.dLocalPixelSizeY);
	}
}

int CResolutionProject::GetNumberOfTargets() const
{
	return m_listOfTargets.size();
}

void CResolutionProject::ClearListOfTargets()
{
	for(size_t i = 0; i < m_listOfTargets.size(); i++)
		delete m_listOfTargets[i];
	
	m_listOfTargets.clear();
}

void CResolutionProject::DeleteTarget(int i)
{
	delete m_listOfTargets[i];
	m_listOfTargets.erase(m_listOfTargets.begin() + i);
}

void CResolutionProject::AddTarget(CTarget* target)
{
	m_listOfTargets.push_back(target);
}

void CResolutionProject::AddEdgeTarget(CEdgeTargetNew *target)
{
	m_listOfTargets.push_back(target);
}

void CResolutionProject::AddLineTarget(CLineTarget *target)
{
	m_listOfTargets.push_back(target);
}

void CResolutionProject::AddLineTarget(const PARE &point)
{
	if(GetDocument())
	{
		Cdfi_Layers * pLayers = GetDocument()->GetLayers();
		CLineTarget *target = new CLineTarget;

		AutoLineDetection(point, target, (int)( min( (m_projectParams.dDefaultTargetWidth / 2) , (m_projectParams.dDefaultTargetLength / 2))));

		target->SetType(LINE_TARGET);

		AddLineTarget(target);
	}
}

void CResolutionProject::AddEdgeTarget(const PARE &point)
{
	if(GetDocument())
	{
		Cdfi_Layers * pLayers = GetDocument()->GetLayers();
		CEdgeTargetNew *target = new CEdgeTargetNew;

		AutoEdgeDetection(point, target, (int)( min( (m_projectParams.dDefaultTargetWidth / 2) , (m_projectParams.dDefaultTargetLength / 2))));

		target->SetType(EDGE_TARGET);

		AddEdgeTarget(target);
	}
}

void CResolutionProject::CalculateDisp(double &dDisp, CTarget *target) const
{
	ESF_PARAMETERS_NEW *ESFparam;
	int nPoints = 0;

	if(target->GetType() == LINE_TARGET)
	{
		nPoints = ((CLineTarget*)target)->GetNumberOfPoints();
		ESFparam = ((CLineTarget*)target)->GetESFparameters();
	}
	else 
	{
		nPoints = ((CEdgeTargetNew*)target)->GetNumberOfPoints();
		ESFparam = ((CEdgeTargetNew*)target)->GetESFparameters();
	}

	double dDisp1 = 0;
	double dMean1 = 0;
	int nDisp1 = 0;
	double dDisp2 = 0;
	int nDisp2 = 0;
	double dMean2 = 0;
	
	double x, value;

	for(int i = 0; i < nPoints; i++)
	{
		INTERPOLATION_POINT point;
		if(target->GetType() == LINE_TARGET)
		{
			x = ((CLineTarget*)target)->GetInterpolationPoint(i)->x;
			value = ((CLineTarget*)target)->GetInterpolationPoint(i)->value;
		}
		else 
		{
			x = ((CEdgeTargetNew*)target)->GetInterpolationPoint(i)->x;
			value = ((CEdgeTargetNew*)target)->GetInterpolationPoint(i)->value;
		}
		
		if(x < ESFparam->x0 - ESFparam->R)
		{
			dDisp1 += value * value;
			dMean1 += value;
			nDisp1++;
		}
		else if(x > ESFparam->x0 + ESFparam->R)
		{
			dDisp2 += value * value;
			dMean2 += value;
			nDisp2++;
		}
	}

	dMean1 /= nDisp1;
	dMean2 /= nDisp2;

	dDisp1 = sqrt((dDisp1 - nDisp1 * dMean1 * dMean1) / (nDisp1 - 1));
	dDisp2 = sqrt((dDisp2 - nDisp2 * dMean2 * dMean2) / (nDisp2 - 1));

	dDisp = max(dDisp1, dDisp2);
}

bool CResolutionProject::SetValidStatus(CTarget *target) const
{
	if(target->GetType() == LINE_TARGET)
	{
		CLineTarget *lineTarget = (CLineTarget *)target;

		ESF_PARAMETERS_NEW *ESFparam = lineTarget->GetESFparameters();

		double dDisp = 0;
		CalculateDisp(dDisp, target);

		target->SetValidStatus(true);

		if(	lineTarget->GetNumberOfPoints() == 0 ||
			ESFparam->R < std::numeric_limits<double>::epsilon() ||
			ESFparam->minValue < 0 || 
			ESFparam->maxValue < 0 ||
			ESFparam->minValue == std::numeric_limits<float>::infinity() ||
			ESFparam->maxValue == std::numeric_limits<float>::infinity() ||
			ESFparam->minValue == -std::numeric_limits<float>::infinity() ||
			ESFparam->maxValue == -std::numeric_limits<float>::infinity() ||
			ESFparam->minValue == std::numeric_limits<double>::quiet_NaN() ||
			ESFparam->maxValue == std::numeric_limits<double>::quiet_NaN() ||
			ESFparam->minValue == -std::numeric_limits<double>::quiet_NaN() ||
			ESFparam->maxValue == -std::numeric_limits<double>::quiet_NaN() ||

			ESFparam->minValue == std::numeric_limits<float>::quiet_NaN() || 
			ESFparam->maxValue == std::numeric_limits<float>::quiet_NaN() || 
			ESFparam->minValue == -std::numeric_limits<float>::quiet_NaN() || 
			ESFparam->maxValue == -std::numeric_limits<float>::quiet_NaN() ||
			abs(ESFparam->x0 - ESFparam->x_min) < (ESFparam->x_max - ESFparam->x_min) / 20 ||
			abs(ESFparam->x0 - ESFparam->x_max) < (ESFparam->x_max - ESFparam->x_min) / 20 
			)
		{
			target->SetValidStatus(false);
			return false;
		}

		else if(ESFparam->maxValue - ESFparam->minValue < (m_projectParams.CSFparamsSNR.dMax - m_projectParams.CSFparamsSNR.dMin) * m_projectParams.dMinTargetRangePercent / 100)
		{
			target->SetValidStatus(false);
			return false;
		}
		else if(dDisp > (ESFparam->maxValue - ESFparam->minValue) * m_projectParams.dDispersePercent/ 100)
		{
			target->SetValidStatus(false);
			return false;
		}
		else if( ESFWidth(ESFparam) > (ESFparam->x_max - ESFparam->x_min) )
		{
			target->SetValidStatus(false);
			return false;
		}
		else
		{
			target->SetValidStatus(true);
			return true;
		}	
	}
	else if(target->GetType() == EDGE_TARGET)
	{
		CEdgeTargetNew *edgeTarget = (CEdgeTargetNew *)target;
	
		ESF_PARAMETERS_NEW *ESFparam = edgeTarget->GetESFparameters();

		double dDisp = 0;
		CalculateDisp(dDisp, target);
		if(	edgeTarget->GetNumberOfPoints() == 0 ||
			ESFparam->R < std::numeric_limits<double>::epsilon() ||
			ESFparam->minValue < 0 || 
			ESFparam->maxValue < 0 ||
			ESFparam->minValue == std::numeric_limits<float>::infinity() ||
			ESFparam->maxValue == std::numeric_limits<float>::infinity() ||
			ESFparam->minValue == -std::numeric_limits<float>::infinity() ||
			ESFparam->maxValue == -std::numeric_limits<float>::infinity() ||
			ESFparam->minValue == std::numeric_limits<double>::quiet_NaN() ||
			ESFparam->maxValue == std::numeric_limits<double>::quiet_NaN() ||
			ESFparam->minValue == -std::numeric_limits<double>::quiet_NaN() ||
			ESFparam->maxValue == -std::numeric_limits<double>::quiet_NaN() ||

			ESFparam->minValue == std::numeric_limits<float>::quiet_NaN() || 
			ESFparam->maxValue == std::numeric_limits<float>::quiet_NaN() || 
			ESFparam->minValue == -std::numeric_limits<float>::quiet_NaN() || 
			ESFparam->maxValue == -std::numeric_limits<float>::quiet_NaN() ||
			abs(ESFparam->x0 - ESFparam->x_min) < (ESFparam->x_max - ESFparam->x_min) / 20 ||
			abs(ESFparam->x0 - ESFparam->x_max) < (ESFparam->x_max - ESFparam->x_min) / 20 
			)
		{
			target->SetValidStatus(false);
			return false;
		}
		else if(ESFparam->maxValue - ESFparam->minValue < (m_projectParams.CSFparamsSNR.dMax - m_projectParams.CSFparamsSNR.dMin) * m_projectParams.dMinTargetRangePercent / 100)
		{
			target->SetValidStatus(false);
			return false;
		}
		else if(dDisp > (ESFparam->maxValue - ESFparam->minValue) * m_projectParams.dDispersePercent/ 100)
		{
			target->SetValidStatus(false);
			return false;
		}
		else if( ESFWidth(ESFparam) > (ESFparam->x_max - ESFparam->x_min) / 1.5)
		{
			target->SetValidStatus(false);
			return false;
		}
		else
		{
			target->SetValidStatus(true);
			return true;
		}	
	}
	return true;
}

//вычисление разрешения по указанному тест-объекту
void CResolutionProject::CalculateResolution(CEdgeTargetNew *target) const
{
	ESF_PARAMETERS_NEW ESFparam;

	if(GetRaster())
	{
		target->SetLocalPixelSizeX(m_projectParams.dLocalPixelSizeX);
		target->SetLocalPixelSizeY(m_projectParams.dLocalPixelSizeY);
		GetESFapproximationParamByTarget(target, ESFparam);
	}
	else
	{
		vector<INTERPOLATION_POINT> points;
		
		int nPoints = target->GetNumberOfPoints();
		for(int i = 0; i < nPoints; i++)
		{
			INTERPOLATION_POINT point;
			point.x = target->GetInterpolationPoint(i)->x;
			point.value = target->GetInterpolationPoint(i)->value;

			points.push_back(point);
			
		}
		if(nPoints > 0)
			GetESFapproximationParamByPoints( points, ESFparam);
	}

	if(target->GetNumberOfPoints() == 0 )
	{
		target->SetPixelResolution(0);
		target->SetMeterResolution(0);
		return;
	}

	SetValidStatus(target);

	double dResolution = 0;
	double dMTFatNyquist = 0;

	int N = 256;

	double dFreqMax = 1.0;

	kiss_fft_cpx * MTFresult;
	MTFresult = new kiss_fft_cpx[N];
	
	CalculateMTF(ESFparam, dFreqMax, N, MTFresult);
	double dMTFvalueMax = sqrt(MTFresult[0].i * MTFresult[0].i + MTFresult[0].r * MTFresult[0].r);

	for(int i = 0; i < N; i ++ )
	{
		double dCSFValue = CSFValue(i * dFreqMax / N);
		
		if(dCSFValue == -1)
		{
			target->SetPixelResolution(0);
			target->SetMeterResolution(0);
			
			dMTFatNyquist = (sqrt(MTFresult[N/2].i * MTFresult[N/2].i + MTFresult[N/2].r * MTFresult[N/2].r) / dMTFvalueMax);
			target->SetMTFatNyquist(dMTFatNyquist);

			delete[] MTFresult;
			return;
		}

		double dMTFvalue = (sqrt(MTFresult[i].i * MTFresult[i].i + MTFresult[i].r * MTFresult[i].r) / dMTFvalueMax); 
		if(dCSFValue > dMTFvalue)
		{
			dResolution = 1. / (i * dFreqMax / N);
			break;
		}
	}

	dMTFatNyquist = (sqrt(MTFresult[N/2].i * MTFresult[N/2].i + MTFresult[N/2].r * MTFresult[N/2].r) / dMTFvalueMax);

	//////////////////////////////////////////
	//////////////////////////////////////////
	
	double dPixelResolution = dResolution / (2);
	target->SetPixelResolution(dPixelResolution);
		
	double dX = dPixelResolution * cos(target->GetAngle()) * m_projectParams.dPixelSizeX;
	double dY = dPixelResolution * sin(target->GetAngle()) * m_projectParams.dPixelSizeY;

	double dMetersResolution = sqrt(dX * dX + dY * dY);
	target->SetMeterResolution(dMetersResolution);

	target->SetMTFatNyquist(dMTFatNyquist);

	delete[] MTFresult;
}

void CResolutionProject::CalculateResolution(CLineTarget *target) const
{
	ESF_PARAMETERS_NEW ESFparam;

	if(GetRaster())
	{
		target->SetLocalPixelSizeX(m_projectParams.dLocalPixelSizeX);
		target->SetLocalPixelSizeY(m_projectParams.dLocalPixelSizeY);
		GetLSFapproximationParamByTarget(target, ESFparam);
	}
	else
	{
		vector<INTERPOLATION_POINT> points;
		for(int i = 0; i < target->GetNumberOfPoints(); i++)
		{
			INTERPOLATION_POINT point;
			point.x = target->GetInterpolationPoint(i)->x;
			point.value = target->GetInterpolationPoint(i)->value;

			points.push_back(point);
		}
		if(target->GetNumberOfPoints() > 0)
			GetLSFapproximationParamByPoints( points, ESFparam, target->IsBrightLine(), TRUE);
	}

	if(target->GetNumberOfPoints() == 0 )
	{
		target->SetPixelResolution(0);
		target->SetMeterResolution(0);
		return;
	}

	SetValidStatus(target);

	double dResolution = 0;

	int N = 256;

	double dFreqMax = 1.0;
	double dMTFatNyquist = 0;

	kiss_fft_cpx * MTFresult;
	MTFresult = new kiss_fft_cpx[N];
	
	CalculateMTF(ESFparam, dFreqMax, N, MTFresult);
	double dMTFvalueMax = sqrt(MTFresult[0].i * MTFresult[0].i + MTFresult[0].r * MTFresult[0].r);

	for(int i = 0; i < N; i ++ )
	{
		double dCSFValue = CSFValue(i * dFreqMax / N);
		if(dCSFValue == -1)
		{
			target->SetPixelResolution(0);
			target->SetMeterResolution(0);
			
			dMTFatNyquist = (sqrt(MTFresult[N/2].i * MTFresult[N/2].i + MTFresult[N/2].r * MTFresult[N/2].r) / dMTFvalueMax);
			target->SetMTFatNyquist(dMTFatNyquist);


			delete[] MTFresult;
			return;
		}
		double dMTFvalue = (sqrt(MTFresult[i].i * MTFresult[i].i + MTFresult[i].r * MTFresult[i].r) / dMTFvalueMax); 
		if(dCSFValue > dMTFvalue)
		{
			dResolution = 1. / (i * dFreqMax / N);
			break;
		}
	}

	dMTFatNyquist = (sqrt(MTFresult[N/2].i * MTFresult[N/2].i + MTFresult[N/2].r * MTFresult[N/2].r) / dMTFvalueMax);

	//////////////////////////////////////////
	//////////////////////////////////////////
	
	double dPixelResolution = dResolution / (2);
	target->SetPixelResolution(dPixelResolution);
		
	double dX = dPixelResolution * cos(target->GetAngle()) * m_projectParams.dPixelSizeX;
	double dY = dPixelResolution * sin(target->GetAngle()) * m_projectParams.dPixelSizeY;

	double dMetersResolution = sqrt(dX * dX + dY * dY);
	target->SetMeterResolution(dMetersResolution);

	target->SetMTFatNyquist(dMTFatNyquist);

	delete[] MTFresult;
}

void CResolutionProject::CalculateResolutionEdge(const PARE& point) const
{
	double resolution = 0;

	CEdgeTargetNew *target = new CEdgeTargetNew;

	AutoEdgeDetection(point, target, 20);

	CalculateResolution( target);
}

void CResolutionProject::FillListOfPoints(CTarget *target, vector<INTERPOLATION_POINT>& resultingListOfPoints) const
{
	////////////////////////////////////////////////////////////////
	// определяем описывающий прямоугольник для тест-объекта
	////////////////////////////////////////////////////////////////

	CIC_PolygonD targetPolygon;
	target->GetPolygon(targetPolygon);

	targetPolygon.UpdateBoundRectGeo();

	CIC_Rect3DD *targetBoundRectGeo = targetPolygon.GetBoundRect();

	targetPolygon.SetExternal(POLYGON_EXTERNAL);

	if(!GetRaster())
	{
		resultingListOfPoints.clear();
		return;
	}

	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;

	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y; 
	double lowerRight_x = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.x; 
	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;

	// Проверяем, что угловые точки попадают внутрь документа  

	PARE point1 = target->GetCornerPoint1();
	PARE point2 = target->GetCornerPoint2();
	PARE point3 = target->GetCornerPoint3();
	PARE point4 = target->GetCornerPoint4();

	if(	point1.x < upperLeft_x || point1.x > lowerRight_x || point1.y > upperLeft_y || point1.y < lowerRight_y ||
		point2.x < upperLeft_x || point2.x > lowerRight_x || point2.y > upperLeft_y || point2.y < lowerRight_y ||
		point3.x < upperLeft_x || point3.x > lowerRight_x || point3.y > upperLeft_y || point3.y < lowerRight_y ||
		point4.x < upperLeft_x || point4.x > lowerRight_x || point4.y > upperLeft_y || point4.y < lowerRight_y )
	{
		resultingListOfPoints.clear();
		if(target->GetType() == LINE_TARGET)
			((CLineTarget*) target)->ClearInterpolationPoints();
		else if(target->GetType() == EDGE_TARGET)
			((CEdgeTargetNew*) target)->ClearInterpolationPoints();
		
		return;
	}

	IC_RECT64 boundRect;

	//вычисляем координаты описывающего прямоугольника
	boundRect.left = (_int64)((targetBoundRectGeo->left - upperLeft_x) / widthPixel);
	boundRect.right = (_int64)((targetBoundRectGeo->right - upperLeft_x) / widthPixel);
	boundRect.top = (_int64)((upperLeft_y - targetBoundRectGeo->top) / heightPixel);
	boundRect.bottom = (_int64)((upperLeft_y - targetBoundRectGeo->bottom) / heightPixel);

	///////////////////////////////////////////////////////////////////
	/// Считываем из растра необходимые пикселы 
	///////////////////////////////////////////////////////////////////

	LPBYTE pBlock;
	INT64 iImageWidth, iImageHeight, wEnd, hEnd;
	int iBlockWidth, iBlockHeight, hb, wb, h, w;
	DWORD dwOffsetInBlock;
	
	GETVALUEPROC pGetValueProc;
	SETVALUEPROC pSetValueProc;
	double dValue,  dMin,  dMax;
	
	CIC_Rect64 rcBlocks;
	CIC_Rect64 rcImage;
	CIC_Rect3DD rcGeo;

	iImageWidth = GetRaster()->GetWidth();
	iImageHeight = GetRaster()->GetHeight();
	iBlockWidth = GetRaster()->GetBlockWidth();
	iBlockHeight = GetRaster()->GetBlockHeight();
	pGetValueProc = CRasterPROCS::PointerGetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));
	pSetValueProc = CRasterPROCS::PointerSetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));

	GetRaster()->GetImageInfo()->GetDefaultRange(GetRaster()->GetPixelType(m_projectParams.iBand), &dMin, &dMax);
		
	rcBlocks.top	= boundRect.top / iBlockHeight;
	rcBlocks.bottom = boundRect.bottom / iBlockHeight;
	rcBlocks.left	= boundRect.left / iBlockWidth;
	rcBlocks.right	= boundRect.right / iBlockWidth;

	Cdfi_RasterStream rasterStream(GetRaster(), 1);
	INT64 iCurrentHeightPosition, iCurrentWidthPosition;
	double dCurrentGeoHeightPosition, dCurrentGeoWidthPosition;
	
	for(hb = (int)rcBlocks.top; hb <= (int)rcBlocks.bottom; hb++)
	{
		for(wb = (int)rcBlocks.left; wb <= (int)rcBlocks.right; wb++) 
		{
			pBlock = (LPBYTE)rasterStream.Get_lpBlock(0, m_projectParams.iBand, 0, wb, hb);
			ASSERT(pBlock);
			
			wEnd = (wb == rcBlocks.right) ? (boundRect.right - wb*iBlockWidth) : (iBlockWidth - 1);
			hEnd = (hb == rcBlocks.bottom) ? (boundRect.bottom - hb*iBlockHeight) : (iBlockHeight - 1);
			
			for(h = 0; h <= hEnd; h++)
			{
				iCurrentHeightPosition = hb * iBlockHeight + h;
				dCurrentGeoHeightPosition = lowerRight_y + (iImageHeight - 1 - iCurrentHeightPosition) * heightPixel;

				for(w = 0, dwOffsetInBlock = (DWORD)(h*iBlockWidth); w <= wEnd; w++, dwOffsetInBlock++)
				{
					iCurrentWidthPosition = wb * iBlockWidth + w;
					dCurrentGeoWidthPosition = upperLeft_x + iCurrentWidthPosition * widthPixel;
										
					// проверка на то, попадает ли точка в прямоугольную область тест-объекта
					if(targetPolygon.PointBelongPolygon(dCurrentGeoWidthPosition, dCurrentGeoHeightPosition) )
					{
						dValue = pGetValueProc((BYTE*)pBlock, dwOffsetInBlock);
						
						//считаем для данной точки расстояние до линии резкого края/ЛПЭ
						double dEdgeDistance = 0;
						if(target->GetType() == LINE_TARGET)
							dEdgeDistance = ((CLineTarget*)target)->GetOrthoEdgeProjection(dCurrentGeoWidthPosition, dCurrentGeoHeightPosition);
						else if(target->GetType() == EDGE_TARGET)
							dEdgeDistance = ((CEdgeTargetNew*)target)->GetOrthoEdgeProjection(dCurrentGeoWidthPosition, dCurrentGeoHeightPosition);
	
						INTERPOLATION_POINT point;
						point.value = dValue;
						point.x = dEdgeDistance;

						resultingListOfPoints.push_back(point);
					}
				}
			}
		}
	}

	if(target->GetType() == EDGE_TARGET)
		SetCorrectSign(resultingListOfPoints);
}

void CResolutionProject::SetCorrectSign(vector<INTERPOLATION_POINT>& resultingListOfPoints) const
{
	size_t i;

	double dSumMinus = 0;
	int nPointMinus = 0;

	double dSumPlus = 0;
	int nPointPlus = 0;
	
	for(i = 0; i < resultingListOfPoints.size(); i++)
	{
		if(resultingListOfPoints[i].x < 0)
		{
			dSumMinus += resultingListOfPoints[i].value;
			nPointMinus++;
		}
		else
		{
			dSumPlus += resultingListOfPoints[i].value;
			nPointPlus++;
		}
	}

	dSumMinus /= nPointMinus;
	dSumPlus /= nPointPlus;
	
	if(dSumMinus > dSumPlus)
		for(i = 0; i < resultingListOfPoints.size(); i++)
			resultingListOfPoints[i].x *= -1;
}

void CResolutionProject::OutputListOfPoints(const vector<INTERPOLATION_POINT> &list, ofstream &resultingFile) const
{
	if (!resultingFile.is_open()) 
		return;

	size_t iSize = list.size();

	for(size_t i = 0; i < iSize; i++)
	{
		resultingFile << list[i].x << " " << list[i].value;

		resultingFile << "\n";
	}
}

double CResolutionProject::GetFunctionalDistance(const vector<INTERPOLATION_POINT>& interpolationPoints,const ESF_PARAMETERS_NEW &params) const
{
	double dDistance = 0;

	double dDelta = min((params.x_max - params.x_min) / 10, 0.5/*GetLocalPixelSize() / 2*/); //шаг для расчета интеграла

	//double dDelta = min((params.x_max - params.x_min) / 10, GetLocalPixelSizeX() / 2); //шаг для расчета интеграла

	for(size_t i = 0; i < interpolationPoints.size(); i++)
	{
		double temp = interpolationPoints[i].value - ESFValue(params, interpolationPoints[i].x, dDelta);
		 
		dDistance += temp * temp;
	}
	return dDistance;
}

double CResolutionProject::GetFunctionalDistanceLSF(const vector<INTERPOLATION_POINT>& interpolationPoints,const ESF_PARAMETERS_NEW &params, bool bBrightLine) const
{
	double dDistance = 0;

	double dDelta = min((params.x_max - params.x_min) / 10, 0.5/*GetLocalPixelSize() / 2*/); //шаг для расчета интеграла

	for(size_t i = 0; i < interpolationPoints.size(); i++)
	{
		double temp = 0;
		if(bBrightLine)
			temp = (interpolationPoints[i].value - params.minValue) / (params.maxValue - params.minValue) - LSFValue(params, interpolationPoints[i].x - params.x0) / LSFValue(params, 0);
		else
			temp = (interpolationPoints[i].value - params.minValue) / (params.maxValue - params.minValue) - (1 - LSFValue(params, interpolationPoints[i].x - params.x0) / LSFValue(params, 0));
		dDistance += temp * temp;
	}

	return dDistance;
}

void CResolutionProject::getMinMaxX(const vector<INTERPOLATION_POINT> &list, double& xMin, double& xMax) const
{
	xMin = list[0].x;
	xMax = list[0].x;

	for(size_t i = 0; i < list.size(); i++)
	{
		xMin = min(xMin, list[i].x);
		xMax = max(xMax, list[i].x);
	}
}

double CResolutionProject::DoStep(	const vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, 
							const double dInitialDistance, double &dResultingDistance, const double dStep, const int iStepIndex) const
{
	double dInitialX = param.x0;
	double dInitialR = param.R;
	double dInitialN = param.N;

	switch(iStepIndex)
	{
		case STEP_X:
			param.x0 += dStep;
		break;
		case STEP_R:
			param.R += dStep;
		break;
		case STEP_N:
			param.N += dStep;
	}

	dResultingDistance = GetFunctionalDistance(points, param);

	if(iStepIndex == STEP_R)
	{
		if(param.R < 0)
		{
			param.x0 = dInitialX;
			param.R = dInitialR;
			param.N = dInitialN;
			dResultingDistance = dInitialDistance;
			
			return 0;
		}
	}

	if(iStepIndex == STEP_N)
	{
		if(param.N < 1.7)
		{
			param.x0 = dInitialX;
			param.R = dInitialR;
			param.N = dInitialN;
			dResultingDistance = dInitialDistance;
			
			return 0;
		}
	}

	if(dResultingDistance < dInitialDistance)
		return dStep;
	else
	{
		switch(iStepIndex)
		{
			case STEP_X:
				param.x0 -= 2 * dStep;
			break;
			case STEP_R:
				param.R -= 2 * dStep;
			break;
			case STEP_N:
				param.N -= 2 * dStep;
		}

		if(iStepIndex == STEP_R)
		{
			if(param.R < 0)
			{
				param.x0 = dInitialX;
				param.R = dInitialR;
				param.N = dInitialN;
				dResultingDistance = dInitialDistance;
			
				return 0;
			}
		}

		if(iStepIndex == STEP_N)
		{
			if(param.N < 1.7)
			{
				param.x0 = dInitialX;
				param.R = dInitialR;
				param.N = dInitialN;
				dResultingDistance = dInitialDistance;
			
				return 0;
			}
		}

		dResultingDistance = GetFunctionalDistance(points, param);

		if(dResultingDistance < dInitialDistance)
			return -dStep;
		else
		{
			param.x0 = dInitialX;
			param.R = dInitialR;
			param.N = dInitialN;
			dResultingDistance = dInitialDistance;
			
			return 0;
		}
	}
	return 0;
}
 
double CResolutionProject::DoStepLSF(	const vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, 
							const double dInitialDistance, double &dResultingDistance, const double dStep, const int iStepIndex, bool bBrightLine) const
{
	double dInitialX = param.x0;
	double dInitialR = param.R;
	double dInitialN = param.N;

	switch(iStepIndex)
	{
		case STEP_X:
			param.x0 += dStep;
		break;
		case STEP_R:
			param.R += dStep;
		break;
		case STEP_N:
			param.N += dStep;
	}

	dResultingDistance = GetFunctionalDistanceLSF(points, param, bBrightLine);

	if(iStepIndex == STEP_R)
	{
		if(param.R < 0)
		{
			param.x0 = dInitialX;
			param.R = dInitialR;
			param.N = dInitialN;
			dResultingDistance = dInitialDistance;
			
			return 0;
		}
	}

	if(iStepIndex == STEP_N)
	{
		if(param.N < 1.7)
		{
			param.x0 = dInitialX;
			param.R = dInitialR;
			param.N = dInitialN;
			dResultingDistance = dInitialDistance;
			
			return 0;
		}
	}

	if(dResultingDistance < dInitialDistance)
		return dStep;
	else
	{
		switch(iStepIndex)
		{
			case STEP_X:
				param.x0 -= 2 * dStep;
			break;
			case STEP_R:
				param.R -= 2 * dStep;
			break;
			case STEP_N:
				param.N -= 2 * dStep;
		}

		if(iStepIndex == STEP_R)
		{
			if(param.R < 0)
			{
				param.x0 = dInitialX;
				param.R = dInitialR;
				param.N = dInitialN;
				dResultingDistance = dInitialDistance;
			
				return 0;
			}
		}

		if(iStepIndex == STEP_N)
		{
			if(param.N < 1.7)
			{
				param.x0 = dInitialX;
				param.R = dInitialR;
				param.N = dInitialN;
				dResultingDistance = dInitialDistance;
			
				return 0;
			}
		}

		dResultingDistance = GetFunctionalDistanceLSF(points, param, bBrightLine);

		if(dResultingDistance < dInitialDistance)
			return -dStep;
		else
		{
			param.x0 = dInitialX;
			param.R = dInitialR;
			param.N = dInitialN;
			dResultingDistance = dInitialDistance;
			
			return 0;
		}
	}
	return 0;
}

void CResolutionProject::GetESFapproximationParamByPoints(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param/*, vector<INTERPOLATION_POINT> &outputPoints*/ ) const
{
	//инициализация начальных параметров
	//строим начальное приближение гауссвовой функцией:
	GetInitialESFParameters(points, param/*, outputPoints*/);
	
	////////////// вывод в файл ///////////////////////////////////////////////////
	
	//ofstream resultingFile;
	//resultingFile.open("C:\\temp\\ESFInitialApproximation.txt", std::ios::out);
	//OutputESFApproximationResults(param, resultingFile);
	//resultingFile.close();
	
	///////////////////////////////////////////////////////////////////////////////

	ASSERT(param.R > 0);

	if(points.size() == 0)
		return;

	double dCurrentStepX = 1.0;
	double dCurrentStepR = 0.25;
	double dCurrentStepN = 0.3;
	
	double dStepX = 0;
	double dStepR = 0;
	double dStepN = 0;
	
	double dInitialDistance = 0;
	double dCurrentDistance = 0;

	int nIterations = 0;

	for(int iNumIterations = 0; iNumIterations < 10 && nIterations < 150/*10*/; ) //while(dStep > 0.007)
	{
		//уточняем значения верхней и нижней полок
		//AdjustMinMaxValue(points, param);

		dInitialDistance = GetFunctionalDistance(points, param);

		//делаем пробные шаги по каждому из параметров
		dStepR = DoStep(points, param, dInitialDistance, dCurrentDistance, dCurrentStepR, STEP_R);
		dStepN = DoStep(points, param, dCurrentDistance, dCurrentDistance, dCurrentStepN, STEP_N);
		dStepX = DoStep(points, param, dCurrentDistance, dCurrentDistance, dCurrentStepX, STEP_X);
		
		if(		dStepX == 0 &&
				dStepN == 0 &&
				dStepR == 0)
		{
			dCurrentStepX /= 2; // в случае, если пробные шаги по всем направлениям не являются успешными - уменьшаем шаг вдвое
			dCurrentStepR /= 2;
			dCurrentStepN /= 2;

			iNumIterations++;
		} 
		else
		{
			dInitialDistance = dCurrentDistance;

			//делаем шаги в том же направлении до тех пор, пока шаг является успешным
			for(bool bSuccessful = true; bSuccessful;)
			{
				//AdjustMinMaxValue(points, param);
				
				nIterations++;

				param.x0 += dStepX;
				param.R += dStepR;
				param.N += dStepN;

				if(param.R < 0 || abs(param.R) < std::numeric_limits<double>::epsilon() || param.N < 1.7)
				{
					param.x0 -= dStepX;
					param.R  -= dStepR;
					param.N  -= dStepN;
					
					bSuccessful = false;
				}
				else
				{

					dCurrentDistance = GetFunctionalDistance(points, param);

					if(dCurrentDistance < dInitialDistance)
						dInitialDistance = dCurrentDistance;
					else
					{
						param.x0 -= dStepX;
						param.R  -= dStepR;
						param.N  -= dStepN;
	
						bSuccessful = false;
					}
				}
			}
		}
	}
	int i = nIterations;
}

void CResolutionProject::GetLSFapproximationParamByPoints(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param, bool bBrightLine, bool bFirstIteration) const
{
	//инициализация начальных параметров
	//строим начальное приближение гауссвовой функцией:
	
	if(bFirstIteration)
		GetInitialLSFParameters(points, param);
	
	//////////////////////////////////////////////////////////////////////////////

	ASSERT(param.R > 0);

	if(points.size() == 0)
		return;

	double dCurrentStepX = 1.0;
	double dCurrentStepR = 0.25;
	double dCurrentStepN = 0.3;
	
	double dStepX = 0;
	double dStepR = 0;
	double dStepN = 0;
	
	double dInitialDistance = 0;
	double dCurrentDistance = 0;

	int nIterations = 0;

	int nMaxIterations;
	
	if(bFirstIteration)
		nMaxIterations = 50;
	else
		nMaxIterations = 150;


	for(int iNumIterations = 0; iNumIterations < 10 && nIterations < nMaxIterations; )
	{
		dInitialDistance = GetFunctionalDistanceLSF(points, param, bBrightLine);

		//делаем пробные шаги по каждому из параметров
		dStepR = DoStepLSF(points, param, dInitialDistance, dCurrentDistance, dCurrentStepR, STEP_R, bBrightLine);
		dStepN = DoStepLSF(points, param, dCurrentDistance, dCurrentDistance, dCurrentStepN, STEP_N, bBrightLine);
		dStepX = DoStepLSF(points, param, dCurrentDistance, dCurrentDistance, dCurrentStepX, STEP_X, bBrightLine);
		
		if(		dStepX == 0 &&
				dStepN == 0 &&
				dStepR == 0)
		{
			dCurrentStepX /= 2; // в случае, если пробные шаги по всем направлениям не являются успешными - уменьшаем шаг вдвое
			dCurrentStepR /= 2;
			dCurrentStepN /= 2;

			iNumIterations++;
		} 
		else
		{
			dInitialDistance = dCurrentDistance;

			//делаем шаги в том же направлении до тех пор, пока шаг является успешным
			for(bool bSuccessful = true; bSuccessful && nIterations < nMaxIterations;)
			{
				nIterations++;

				param.x0 += dStepX;
				param.R += dStepR;
				param.N += dStepN;

				if(param.R < 0 || abs(param.R) < std::numeric_limits<double>::epsilon() || param.N < 1.4)
				{
					param.x0 -= dStepX;
					param.R  -= dStepR;
					param.N  -= dStepN;
					
					bSuccessful = false;
				}
				else
				{
					dCurrentDistance = GetFunctionalDistanceLSF(points, param, bBrightLine);

					if(dCurrentDistance < dInitialDistance)
						dInitialDistance = dCurrentDistance;
					else
					{
						param.x0 -= dStepX;
						param.R  -= dStepR;
						param.N  -= dStepN;
	
						bSuccessful = false;
					}
				}
			}
		}
	}

	if(bFirstIteration)
	{
		double x1 = param.x0 - 1.2 * param.R;
		double x2 = param.x0 + 1.2 * param.R;

		int nPoints = 0;
		double dMean = 0;

		for(size_t i = 0; i < points.size(); i++)
		{
			if(points[i].x < x1 || points[i].x > x2)
			{
				nPoints++;
				dMean += points[i].value; 
			}
		}

		if(nPoints != 0)
		{
			dMean /= nPoints;
		
			if(bBrightLine)
				param.minValue = dMean;
			else
				param.maxValue = dMean;
			
			GetLSFapproximationParamByPoints(points, param, bBrightLine, FALSE);
		}
	}  
}

void CResolutionProject::OutputESFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const
{
	if (!resultingFile.is_open()) 
		return;
		
	double dStep = 0.1;

	for(double x = param.x_min; x < param.x_max; x += dStep)
	{
		resultingFile << x << " " << ESFValue(param, x, 0.01);

		resultingFile << "\n";
	}
}

void CResolutionProject::OutputLSFPoints(const ESF_PARAMETERS_NEW& param, const vector<INTERPOLATION_POINT> &list, ofstream &resultingFile) const
{
	if (!resultingFile.is_open()) 
		return;

	size_t iSize = list.size();

	for(size_t i = 0; i < iSize; i++)
	{
		resultingFile << list[i].x << " " << (list[i].value - param.minValue) / (param.maxValue - param.minValue) ;

		resultingFile << "\n";
	}
}

void CResolutionProject::GetESFapproximationParamByTarget(CEdgeTargetNew *target, ESF_PARAMETERS_NEW& param ) const
{
	vector<INTERPOLATION_POINT> list;
	FillListOfPoints(target, list);
	
	if(list.size() == 0)
		return;
	
	InterpolatePoints(list);

	GetESFapproximationParamByPoints(list, param);
	
	target->SetInterpolationPoints(list);
	target->SetESFparameters(param);
}

void CResolutionProject::GetLSFapproximationParamByTarget(CLineTarget *target, ESF_PARAMETERS_NEW& param ) const
{
	vector<INTERPOLATION_POINT> list;
	FillListOfPoints(target, list);
	
	if(list.size() == 0)
		return;
	
	InterpolatePointsLSF(list);

	target->SetInterpolationPoints(list);
	target->SetBrightnessType();

	GetLSFapproximationParamByPoints(list, param, target->IsBrightLine(), TRUE);

	target->SetESFparameters(param);
}

double CResolutionProject::ESFValue(const ESF_PARAMETERS_NEW& param, double x, double dDelta)
{
	//  формула аппроксимации 
	//  (ESF - edge spread function; переходная функция):
	//
	//                                                  x
	//												     _                
	//												    |     /          | z - x0 | N  \
	// ESF(x) = minValue + (maxValue - minValue) * S *  | exp | -- 4.6 * | ------ |    |  dz
	//												   _|     \          |   R    |    /
	//           
	//                                                 x_min
	//
	// магический смысл параметра 4.6 - при введение этого параметра R примерно равняется ширине LSF на уровне 0.1 от ее высоты
	//
	
	//  считаем интеграл 
	//              x_max
	//	              _                
	//	             |     /          | z - x0 | N  \                      /          | z - x0 | N \
	//  dIntegral =  | exp | -- 4.6 * | ------ |    |  dz ,    функцию exp | -- 4.6 * | ------ |   |   обозначим  Exp(const ESF_PARAMETERS& param, double z)
	//	            _|     \          |   R    |    /                      \          |   R    |   /
	//          
	//              x_min

	//double dDelta = min((param.x_max - param.x_min) / 10, 2 * m_dPixelSize); //шаг для расчета интеграла
	
	double dIntegral = 0;

	for(double z = param.x_min; z < param.x_max; z += dDelta)
		dIntegral += dDelta * Exp(param, z);
	
	// считаем выражение
	//
	//                x
	//	               _                
	//	              |     /          | z - x0 | N  \          
	//  dXintegral =  | exp | -- 4.6 * | ------ |    |  dz ,    
	//	             _|     \          |   R    |    /          
	//          
	//              x_min

	double dXintegral = 0;

	for(double z = param.x_min; z < x; z += dDelta)
		dXintegral += dDelta * Exp(param, z);

	double result = param.minValue + (param.maxValue - param.minValue) * (dXintegral / dIntegral);

	return result;
}

double CResolutionProject::ESFWidth(const ESF_PARAMETERS_NEW& param)
{
	double dWidth = param.x_max - param.x_min;
	double dStep = param.maxValue - param.minValue;

	for(double x = param.x_min; x < param.x_max; x += dWidth / 20)
	{
		if(ESFValue(param, x, (param.x_max - param.x_min) / 10) > param.minValue + dStep / 100)
			return 2 * (param.x0 - x + dWidth / 20);
	}

	return 0;
}

double CResolutionProject::Exp(const ESF_PARAMETERS_NEW& param, double z)
{
	double result = 0;

	//	           /          | z - x0 | N  \          
	//  Exp =  exp | -- 4.6 * | ------ |    |  
	//	           \          |   R    |    /          
	          
	result = fabsl((z - param.x0) / param.R);
	result = - 4.6 * pow(result, param.N);
	result = exp(result);

	return result;
}
 
double CResolutionProject::Exp0(const ESF_PARAMETERS_NEW& param, double z)// const
{
	double result = 0;

	//	            /          |   z   | N  \          
	//  Exp0 =  exp | -- 4.6 * | ----- |    |  
	//	            \          |   R   |    /          
	          
	result = fabsl(z / param.R);
	result = - 4.6 * pow(result, param.N);
	result = exp(result);

	return result;
}

double CResolutionProject::LSFValue(const ESF_PARAMETERS_NEW& param, double x)// const
{
	//									        /          |   x   | N  \
	// LSF(x) = (maxValue - minValue) * S * exp | -- 4.6 * | ----- |    |  
	//										    \          |   R   |    /
	
	// LSF(x) возвращает функцию, симметричную относительно оси (x = 0) - параметр сдвига x0 не учитывается

	double dDelta = min((param.x_max - param.x_min) / 10, 0.2 /*GetLocalPixelSizeX() / 5*/);
	//dDelta = 0.01; //шаг для расчета интеграла
	double dIntegral = 0; //dIntegral = 1/S

	for(double z = param.x_min; z < param.x_max; z += dDelta)
		dIntegral += dDelta * Exp0(param, z);

	double result = (param.maxValue - param.minValue) * (Exp0(param, x) / dIntegral);

	return result;
}    

void  CResolutionProject::OutputLSFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const
{
	if (!resultingFile.is_open()) 
		return;
		
	double dStep = 0.1;

	for(double x = param.x_min; x < param.x_max; x += dStep)
	{
		resultingFile << x << " " << LSFValue(param, x - param.x0) / LSFValue(param, 0);
		resultingFile << "\n";
	}
}

void CResolutionProject::CalculateMTF(ESF_PARAMETERS_NEW *param, int nESFplots, double dFreqMax, int N, kiss_fft_cpx *out) const
{
	// Для расчета FFT используется библиотека "kissFFT" (http://sourceforge.net/projects/kissfft/)
	
	// ESF_PARAMETERS& param - параметры аналитической аппроксимации LSF
	// kiss_fft_cpx *out - массив для вывода результата (значений амплитуд, т.е. модулей комплексных чисел, получающихся в результате работы FFT)
	// N количество точек в выходном массиве out 
	
	// т.е. частота дискретизации равна dStepFreq = dFreqMax / N (***)

	// (частоты соответствующие выходным значениям равны {0, dFreqMax / N, 2 * dFreqMax / N, ... , dFreqMax} )

	// ( реальный массив, получающийся на выходе алгоритма FFT имеет размер 2N, 
	//   но при вещественных исходных данных последний N точек в нем симметрично повторяют первые N, 
	//   поэтому в массив out передаются только первые N точек )

	// Таким образом, размер входного массива данных должен быть равен 2N
	// Обозначим dStepX интервал между точками LSF (т.е шаг, с которым проходит дискретизация по расстоянию функции LSF)
	// Шаг результирующих значений по частоте (частота дискретизации) равен (http://habrahabr.ru/post/196374/):
	
	//                                               _
	//         dStepFreq = 1 / (2N * dStepX)          \
	//                                                 |
	// С другой стороны, как отмечалось ранее (***)    }   =>   1 / (2N * dStepX) =  dFreqMax / N  =>   dStepX =  1 / (2 * dFreqMax)  (****)
	//                                                 |
	//		   dStepFreq = dFreqMax / N              _/
	
	double dStepFreq = dFreqMax / N;
	double dStepX = 1 / (2 * dFreqMax); //см. (****)
	double dXmin = - N * dStepX;
	
	///  расчет входных данных для алгоритма FFT //
	kiss_fft_cpx *in, *outTmp;
	in = new kiss_fft_cpx[2 * N];
	outTmp = new kiss_fft_cpx[2 * N];
	
	for(int iX = 0; iX < 2 * N; iX++)
	{
		in[iX].r = 0;
		for(int j = 0; j < nESFplots; j++)
			in[iX].r += (float)LSFValue(param[j], dXmin + iX * dStepX);	
		
		in[iX].r /= nESFplots;
		
		in[iX].i = 0;
		outTmp[iX].r = 0;
		outTmp[iX].i = 0;
	}
	//////////////////////////////////////////////

	/////////////  Расчет FFT  ///////////////////
	kiss_fft_cfg cfg;

	if ((cfg = kiss_fft_alloc(2 * N, 0/*is_inverse_fft*/, NULL, NULL)) != NULL)
	{
		kiss_fft(cfg, in, outTmp);
		free(cfg);
	}
	//////////////////////////////////////////////

	for(int iX = 0; iX < N; iX++)
	{
		out[iX].r = outTmp[iX].r;
		out[iX].i = outTmp[iX].i;
	}

	delete[] in;
	delete[] outTmp;
}

void CResolutionProject::CalculateMTF(const ESF_PARAMETERS_NEW& param, double dFreqMax, int N, kiss_fft_cpx *out) const
{
	// Для расчета FFT используется библиотека "kissFFT" (http://sourceforge.net/projects/kissfft/)
	
	// ESF_PARAMETERS& param - параметры аналитической аппроксимации LSF
	// kiss_fft_cpx *out - массив для вывода результата (значений амплитуд, т.е. модулей комплексных чисел, получающихся в результате работы FFT)
	// N количество точек в выходном массиве out 
	
	// т.е. частота дискретизации равна dStepFreq = dFreqMax / N (***)

	// (частоты соответствующие выходным значениям равны {0, dFreqMax / N, 2 * dFreqMax / N, ... , dFreqMax} )

	// ( реальный массив, получающийся на выходе алгоритма FFT имеет размер 2N, 
	//   но при вещественных исходных данных последний N точек в нем симметрично повторяют первые N, 
	//   поэтому в массив out передаются только первые N точек )

	// Таким образом, размер входного массива данных должен быть равен 2N
	// Обозначим dStepX интервал между точками LSF (т.е шаг, с которым проходит дискретизация по расстоянию функции LSF)
	// Шаг результирующих значений по частоте (частота дискретизации) равен (http://habrahabr.ru/post/196374/):
	
	//                                               _
	//         dStepFreq = 1 / (2N * dStepX)          \
	//                                                 |
	// С другой стороны, как отмечалось ранее (***)    }   =>   1 / (2N * dStepX) =  dFreqMax / N  =>   dStepX =  1 / (2 * dFreqMax)  (****)
	//                                                 |
	//		   dStepFreq = dFreqMax / N              _/
	
	double dStepFreq = dFreqMax / N;
	double dStepX = 1 / (2 * dFreqMax); //см. (****)
	double dXmin = - N * dStepX;
	
	///  расчет входных данных для алгоритма FFT //
	kiss_fft_cpx *in, *outTmp;
	in = new kiss_fft_cpx[2 * N];
	outTmp = new kiss_fft_cpx[2 * N];
	
	for(int iX = 0; iX < 2 * N; iX++)
	{
		in[iX].r = (float)LSFValue(param, dXmin + iX * dStepX);
		in[iX].i = 0;
		outTmp[iX].r = 0;
		outTmp[iX].i = 0;
	}
	//////////////////////////////////////////////

	/////////////  Расчет FFT  ///////////////////
	kiss_fft_cfg cfg;

	if ((cfg = kiss_fft_alloc(2 * N, 0, NULL, NULL)) != NULL)
	{
		kiss_fft(cfg, in, outTmp);
		free(cfg);
	}
	//////////////////////////////////////////////

	for(int iX = 0; iX < N; iX++)
	{
		out[iX].r = outTmp[iX].r;
		out[iX].i = outTmp[iX].i;
	}

	delete[] in;
	delete[] outTmp;
}

void  CResolutionProject::OutputMTFApproximationResults(const ESF_PARAMETERS_NEW& param, ofstream &resultingFile) const
{
	if (!resultingFile.is_open()) 
		return;
	
	int N = 128/*256*/;
	double dFreqMax = 3. / GetLocalPixelSizeX();
	
	kiss_fft_cpx *out;
	out = new kiss_fft_cpx[N];
	
	CalculateMTF(param, dFreqMax, N, out);

	double dFreqStep = dFreqMax / N; // частота дискретизации

	double dMaxMTF = sqrt (out[0].r * out[0].r + out[0].i * out[0].i);

	for(int iFreq = 0; iFreq < N; iFreq++)
	{
		double r = out[iFreq].r;
		double i = out[iFreq].i;
		resultingFile << iFreq * dFreqStep  << " " << sqrt(r * r + i * i) / dMaxMTF ;

		resultingFile << "\n";
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// TEST CODE FOR DIRECT MTF CALCULATION WITHOUT FFT ////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/*resultingFile << "\n\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";

	double dMaxMTF = MTFValue(param, 0).Abs();

	double dStep = dFreqStep;

	for(double dFrequency = 0; dFrequency < 3. / m_dPixelSize; dFrequency += dStep)
	{
		resultingFile << dFrequency << " " << MTFValue(param, dFrequency).Abs() / dMaxMTF;

		resultingFile << "\n";
	}*/

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 		
	delete[] out;
}

void CResolutionProject::OutputCSFApproximationResults(ofstream &resultingFile) const
{
	if (!resultingFile.is_open()) 
		return;
		
	double dStep = 0.2;

	for(double dFrequency = 0; dFrequency < 2 / GetLocalPixelSizeX(); dFrequency += dStep)
	{
		resultingFile << dFrequency << " " << CSFValue(dFrequency);
		resultingFile << "\n";
	}
}

void CResolutionProject::AdjustMinMaxValue(const vector<INTERPOLATION_POINT>& list, ESF_PARAMETERS_NEW& param) const
{
	// находим значения верхней и нижней полок по методу наименьших квадратов
	// представляем ERF(x) = maxValue + minValue * g(x), где
	
	//                 x
	//				    _                
	//			       |     /          | z - x0 | N  \
	//     g(x) = S *  | exp | -- 4.6 * | ------ |    |  dz  ,
	//				  _|     \          |   R    |    /
	//           
	//               x_min

	double sum_ERF = 0;
	double sum_g = 0;
	double sum_g2 = 0;
	double sum_ERF_g = 0;

	int N = list.size(); //количество точек аппроксимации в массиве list
	
	if(N == 0)
	{
		param.minValue = 0;
		param.maxValue = 0;
		return;
	}

	for(int  i = 0; i < N; i++)
	{
		sum_ERF += list[i].value;
		double tmp = g(param, list[i].x);
		sum_g += tmp;
		sum_g2 += tmp * tmp;
		sum_ERF_g += list[i].value * tmp;
	}

	//формула аналогичная формуле для нахождения коэффициентов линейной регресси по методу наименьших квадратов

	double dDeltaH = (N * sum_ERF_g - sum_ERF * sum_g) / (N * sum_g2 - sum_g * sum_g);
	 
	param.minValue = (sum_ERF - dDeltaH * sum_g) / N;

	param.maxValue = param.minValue + dDeltaH;
}
 
double CResolutionProject::g(const ESF_PARAMETERS_NEW& param, double x) const
{
	//                                                  x
	//												     _                
	//												    |     /    | z - x0 | N  \
	// ESF(x) = minValue + (maxValue - minValue) * S *  | exp | -- | ------ |    |  dz
	//												   _|     \    |   R    |    /
	//           
	//                                                 x_min
	//                 x
	//				    _                
	//			       |     /    | z - x0 | N  \
	//     g(x) = S *  | exp | -- | ------ |    |  dz  ,
	//				  _|     \    |   R    |    /
	//           
	//               x_min

	double dDelta = min((param.x_max - param.x_min) / 10, 0.2 /*GetLocalPixelSizeX() / 5*/); //шаг для расчета интеграла

	return (ESFValue(param, x, dDelta) - param.minValue) / (param.maxValue - param.minValue);
}

CIC_ComplexDouble CResolutionProject::MTFValue(const ESF_PARAMETERS_NEW& param, double x) const
{
	//вычисляем преобразование Фурье от LSF
	double dStep = 0.1;

	CIC_ComplexDouble result;
	result.im = 0;
	result.re = 0;
	
	for(double x_cur = - (param.x_max - param.x_min) / 2 ; x_cur < (param.x_max - param.x_min) / 2; x_cur += dStep)
	{
		result.re += cos(2 * M_PI * x * x_cur) * LSFValue(param, x_cur) * dStep;
		result.im += - sin(2 * M_PI * x * x_cur) * LSFValue(param, x_cur) * dStep;
	}
	return result;
}

double CResolutionProject::ApproximateAngleDetection(const PARE &point, double Rpix) const
{
	double dResultAngle = 0;
	double dMaxDelta = 0;
	double dCurrentDelta = 0;
	
	double dSum1 = 0;
	double dSum2 = 0;

	CVector tmpVector;
		
	LPBYTE pBlock;
	INT64 iImageWidth, iImageHeight, wEnd, hEnd;
	int iBlockWidth, iBlockHeight, hb, wb, h, w;
	DWORD dwOffsetInBlock;
	
	if(!GetRaster())
		return -1;

	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;

	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y; 
	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;

	//DWORD dwPixelSize, dwBlockSize, dwOffsetInBlock;
	GETVALUEPROC pGetValueProc;
	SETVALUEPROC pSetValueProc;
	double dValue,  dMin,  dMax;
	
	CIC_Rect64 rcBlocks;
	CIC_Rect64 rcImage;
	CIC_Rect3DD rcGeo;

	iImageWidth = GetRaster()->GetWidth();
	iImageHeight = GetRaster()->GetHeight();
	iBlockWidth = GetRaster()->GetBlockWidth();
	iBlockHeight = GetRaster()->GetBlockHeight();

	Cdfi_Raster *pRaster = GetRaster();
	int iBand = m_projectParams.iBand;

	pGetValueProc = CRasterPROCS::PointerGetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));
	pSetValueProc = CRasterPROCS::PointerSetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));

	GetRaster()->GetImageInfo()->GetDefaultRange(GetRaster()->GetPixelType(m_projectParams.iBand), &dMin, &dMax);
		
	IC_RECT64 boundRect;

	//вычисляем координаты описывающего прямоугольника
	boundRect.left   = (_int64)((point.x - upperLeft_x) / widthPixel - Rpix);
	boundRect.right  = (_int64)((point.x - upperLeft_x) / widthPixel + Rpix);
	boundRect.top    = (_int64)((upperLeft_y - point.y) / heightPixel - Rpix);
	boundRect.bottom = (_int64)((upperLeft_y - point.y) / heightPixel + Rpix);

	rcBlocks.top	= boundRect.top / iBlockHeight;
	rcBlocks.bottom = boundRect.bottom / iBlockHeight;
	rcBlocks.left	= boundRect.left / iBlockWidth;
	rcBlocks.right	= boundRect.right / iBlockWidth;

	Cdfi_RasterStream rasterStream(GetRaster(), 1);
	INT64 iCurrentHeightPosition, iCurrentWidthPosition;
	double dCurrentGeoHeightPosition, dCurrentGeoWidthPosition;
	
	for(double angle = 0; angle < M_PI; angle += M_PI / 6)
	{
		dSum1 = 0;
		dSum2 = 0;

		for(hb = (int)rcBlocks.top; hb <= (int)rcBlocks.bottom; hb++)
		{
			for(wb = (int)rcBlocks.left; wb <= (int)rcBlocks.right; wb++) 
			{
				pBlock = (LPBYTE)rasterStream.Get_lpBlock(0, m_projectParams.iBand, 0, wb, hb);
				ASSERT(pBlock);
			
				wEnd = (wb == rcBlocks.right) ? (boundRect.right - wb * iBlockWidth) : (iBlockWidth - 1);
				hEnd = (hb == rcBlocks.bottom) ? (boundRect.bottom - hb * iBlockHeight) : (iBlockHeight - 1);
			
				for(h = 0; h <= hEnd; h++)
				{
					iCurrentHeightPosition = hb * iBlockHeight + h;
					dCurrentGeoHeightPosition = lowerRight_y + (iImageHeight - 1 - iCurrentHeightPosition) * heightPixel;

					for(w = 0, dwOffsetInBlock = (DWORD)(h*iBlockWidth); w <= wEnd; w++, dwOffsetInBlock++)
					{
						iCurrentWidthPosition = wb * iBlockWidth + w;
						dCurrentGeoWidthPosition = upperLeft_x + iCurrentWidthPosition * widthPixel;
								
						double dXPix = (dCurrentGeoWidthPosition - point.x) / widthPixel;
						double dYPix = (dCurrentGeoHeightPosition - point.y) / heightPixel;
						// проверка на то, попадает ли точка в прямоугольную область тест-объекта
						if(sqrt(dXPix * dXPix + dYPix * dYPix) < Rpix)
							//sqrt((dCurrentGeoWidthPosition - point.x) * (dCurrentGeoWidthPosition - point.x) + (dCurrentGeoHeightPosition - point.y) * (dCurrentGeoHeightPosition - point.y)) < R)
						{
							dValue = pGetValueProc((BYTE*)pBlock, dwOffsetInBlock);
							tmpVector.InitializeXY(dCurrentGeoWidthPosition - point.x, dCurrentGeoHeightPosition - point.y);
							double tmpAngle = tmpVector.GetAngle();
							if(tmpAngle > angle || tmpAngle < angle - M_PI) 
								dSum1 += dValue;
							else
								dSum2 += dValue;
						}
					}
				}
				GetRaster()->ReleaseLPBlock(m_projectParams.iBand, 0, wb, hb);
			}
		}

		dCurrentDelta = fabsl(dSum1 - dSum2);
		if(dCurrentDelta > dMaxDelta)
		{
			dMaxDelta = dCurrentDelta;
			dResultAngle = angle;
		}
	}
	//цикл перебора углов с шагом в 30 градусов

	return dResultAngle;
}


double CResolutionProject::ApproximateAngleDetectionLSF(const PARE &point, double Rpix) const
{
	double dResultAngle = 0;
	double dMaxDelta = 0;
	double dCurrentDelta = 0;
	
	double dSum1 = 0;
	double dSum2 = 0;

	CVector tmpVector, tmpVector2;
		
	LPBYTE pBlock;
	INT64 iImageWidth, iImageHeight, wEnd, hEnd;
	int iBlockWidth, iBlockHeight, hb, wb, h, w;
	DWORD dwOffsetInBlock;
	
	if(!GetRaster())
		return -1;

	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;

	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y; 
	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;

	GETVALUEPROC pGetValueProc;
	SETVALUEPROC pSetValueProc;
	double dValue,  dMin,  dMax;
	
	CIC_Rect64 rcBlocks;
	CIC_Rect64 rcImage;
	CIC_Rect3DD rcGeo;

	iImageWidth = GetRaster()->GetWidth();
	iImageHeight = GetRaster()->GetHeight();
	iBlockWidth = GetRaster()->GetBlockWidth();
	iBlockHeight = GetRaster()->GetBlockHeight();

	Cdfi_Raster *pRaster = GetRaster();
	int iBand = m_projectParams.iBand;

	pGetValueProc = CRasterPROCS::PointerGetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));
	pSetValueProc = CRasterPROCS::PointerSetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));

	GetRaster()->GetImageInfo()->GetDefaultRange(GetRaster()->GetPixelType(m_projectParams.iBand), &dMin, &dMax);
		
	IC_RECT64 boundRect;

	//вычисляем координаты описывающего прямоугольника
	boundRect.left   = (_int64)((point.x - upperLeft_x) / widthPixel - Rpix);
	boundRect.right  = (_int64)((point.x - upperLeft_x) / widthPixel + Rpix);
	boundRect.top    = (_int64)((upperLeft_y - point.y) / heightPixel - Rpix);
	boundRect.bottom = (_int64)((upperLeft_y - point.y) / heightPixel + Rpix);

	rcBlocks.top	= boundRect.top / iBlockHeight;
	rcBlocks.bottom = boundRect.bottom / iBlockHeight;
	rcBlocks.left	= boundRect.left / iBlockWidth;
	rcBlocks.right	= boundRect.right / iBlockWidth;

	Cdfi_RasterStream rasterStream(GetRaster(), 1);
	INT64 iCurrentHeightPosition, iCurrentWidthPosition;
	double dCurrentGeoHeightPosition, dCurrentGeoWidthPosition;
	
	bool bBegin = true;
	double dMinDisp = 0;

	//////////////////////////////////////////
	// Определяем темная линия или светлая
	//////////////////////////////////////////
	double dMeanSmall = 0;
	int nPointsSmall = 0;
	double dMeanBig = 0;
	int nPointsBig = 0;
	double dMaxValue = 0;

	for(hb = (int)rcBlocks.top; hb <= (int)rcBlocks.bottom; hb++)
	{
		for(wb = (int)rcBlocks.left; wb <= (int)rcBlocks.right; wb++) 
		{
			pBlock = (LPBYTE)rasterStream.Get_lpBlock(0, m_projectParams.iBand, 0, wb, hb);
			ASSERT(pBlock);
			
			wEnd = (wb == rcBlocks.right) ? (boundRect.right - wb * iBlockWidth) : (iBlockWidth - 1);
			hEnd = (hb == rcBlocks.bottom) ? (boundRect.bottom - hb * iBlockHeight) : (iBlockHeight - 1);
			
			for(h = 0; h <= hEnd; h++)
			{
				iCurrentHeightPosition = hb * iBlockHeight + h;
				dCurrentGeoHeightPosition = lowerRight_y + (iImageHeight - 1 - iCurrentHeightPosition) * heightPixel;

				for(w = 0, dwOffsetInBlock = (DWORD)(h*iBlockWidth); w <= wEnd; w++, dwOffsetInBlock++)
				{
					iCurrentWidthPosition = wb * iBlockWidth + w;
					dCurrentGeoWidthPosition = upperLeft_x + iCurrentWidthPosition * widthPixel;
								
					double dXPix = (dCurrentGeoWidthPosition - point.x) / widthPixel;
					double dYPix = (dCurrentGeoHeightPosition - point.y) / heightPixel;
					// проверка на то, попадает ли точка в прямоугольную область тест-объекта
					if(sqrt(dXPix * dXPix + dYPix * dYPix) < Rpix)
					{
						dValue = pGetValueProc((BYTE*)pBlock, dwOffsetInBlock);
						dMaxValue = max(dValue, dMaxValue);

						if(sqrt(dXPix * dXPix + dYPix * dYPix) < min(Rpix, 5))
						{
							dValue = pGetValueProc((BYTE*)pBlock, dwOffsetInBlock);
							
							dMeanBig += dValue;
							nPointsBig++;

							if(sqrt(dXPix * dXPix + dYPix * dYPix) < 2)
							{
								dMeanSmall += dValue;	
								nPointsSmall++;
							}
						}
					}
				}
			}
			GetRaster()->ReleaseLPBlock(m_projectParams.iBand, 0, wb, hb);
		}
	}

	dMeanSmall /= nPointsSmall;
	dMeanBig /= nPointsBig;

	bool bInverse = false;

	if(dMeanSmall > dMeanBig)
		bInverse = true;
		
	
	for(double angle = 0; angle < M_PI; angle += M_PI / 90)
	{
		dSum1 = 0;
		dSum2 = 0;

		double dSumValues = 0;
		double dSumValues2 = 0;
		int nValues = 0;
		
		for(hb = (int)rcBlocks.top; hb <= (int)rcBlocks.bottom; hb++)
		{
			for(wb = (int)rcBlocks.left; wb <= (int)rcBlocks.right; wb++) 
			{
				pBlock = (LPBYTE)rasterStream.Get_lpBlock(0, m_projectParams.iBand, 0, wb, hb);
				ASSERT(pBlock);
			
				wEnd = (wb == rcBlocks.right) ? (boundRect.right - wb * iBlockWidth) : (iBlockWidth - 1);
				hEnd = (hb == rcBlocks.bottom) ? (boundRect.bottom - hb * iBlockHeight) : (iBlockHeight - 1);
			
				for(h = 0; h <= hEnd; h++)
				{
					iCurrentHeightPosition = hb * iBlockHeight + h;
					dCurrentGeoHeightPosition = lowerRight_y + (iImageHeight - 1 - iCurrentHeightPosition) * heightPixel;

					for(w = 0, dwOffsetInBlock = (DWORD)(h*iBlockWidth); w <= wEnd; w++, dwOffsetInBlock++)
					{
						iCurrentWidthPosition = wb * iBlockWidth + w;
						dCurrentGeoWidthPosition = upperLeft_x + iCurrentWidthPosition * widthPixel;
								
						double dXPix = (dCurrentGeoWidthPosition - point.x) / widthPixel;
						double dYPix = (dCurrentGeoHeightPosition - point.y) / heightPixel;
						// проверка на то, попадает ли точка в прямоугольную область тест-объекта
						if(sqrt(dXPix * dXPix + dYPix * dYPix) < Rpix)
							//sqrt((dCurrentGeoWidthPosition - point.x) * (dCurrentGeoWidthPosition - point.x) + (dCurrentGeoHeightPosition - point.y) * (dCurrentGeoHeightPosition - point.y)) < R)
						{
							dValue = pGetValueProc((BYTE*)pBlock, dwOffsetInBlock);
							tmpVector.InitializeXY(dCurrentGeoWidthPosition - point.x, dCurrentGeoHeightPosition - point.y);

							tmpVector.InitializeXY(dXPix, dYPix);
							tmpVector2.InitializeAngle(angle - M_PI / 2, 1);

							double length = Scalar(tmpVector, tmpVector2);

							if(abs(length) < 2)
							{
								if(!bInverse)
								{
									dSumValues += dValue;
									dSumValues2 += dValue * dValue;
								}
								else
								{
									dSumValues += dMaxValue - dValue;
									dSumValues2 += (dMaxValue - dValue) * (dMaxValue - dValue);
								}
								nValues++;
							}
						}
					}
				}
				GetRaster()->ReleaseLPBlock(m_projectParams.iBand, 0, wb, hb);
			}
		}

		double disp = (dSumValues2 - (dSumValues / nValues) * (dSumValues / nValues)) / (nValues - 1);

		if(disp < dMinDisp || bBegin)
		{
			dMinDisp = disp;
			dResultAngle = angle;
			bBegin = false;
		}
	}
	//цикл перебора углов с шагом в 30 градусов
	//

	return dResultAngle;
}

void CResolutionProject::AutoLineDetection(const PARE &point, CLineTarget* target, int iSearchRadius /*радиус*/) const
{
	//length - длина в пикселах 
	//получаем приблизительное (с точностью в 30 градусов положение резкого края) в промежутке [0, PI]
	if(!GetRaster())
		return;
	
	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;

	double dRadius = widthPixel * iSearchRadius;
		
	CIC_Rect3DD geoRect;
	geoRect.top = point.y + iSearchRadius * heightPixel;
	geoRect.bottom = point.y - iSearchRadius * heightPixel;
	geoRect.left = point.x - iSearchRadius * widthPixel;
	geoRect.right = point.x + iSearchRadius * widthPixel;
	
	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;
	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y;

	INT64 iImageHeight = GetRaster()->GetHeight();

	IC_RECT64 boundRect;

	boundRect.left = (_int64)((point.x  - upperLeft_x) / widthPixel - iSearchRadius);
	boundRect.right = (_int64)((point.x  - upperLeft_x) / widthPixel + iSearchRadius);
	boundRect.top = (_int64)((upperLeft_y - point.y) / heightPixel - iSearchRadius);
	boundRect.bottom = (_int64)((upperLeft_y - point.y) / heightPixel + iSearchRadius);

	//получаем приближенное значение угла

	double angle = ApproximateAngleDetectionLSF(point, (double)iSearchRadius);

	CVector vector;
	vector.InitializeAngle(angle, m_projectParams.dDefaultTargetLength / 2);
	double length = vector.GetLength();

	PARE edgePoint1, edgePoint2;
	edgePoint1.x = point.x + vector.GetX() * widthPixel;
	edgePoint1.y = point.y + vector.GetY() * heightPixel;

	edgePoint2.x = point.x - vector.GetX() * widthPixel;
	edgePoint2.y = point.y - vector.GetY() * heightPixel;

	target->Initialize(edgePoint1, edgePoint2, m_projectParams.dDefaultTargetWidth / 2, m_projectParams.dDefaultTargetWidth / 2, widthPixel, heightPixel);
	
	#ifdef _DEBUG
		//PointsVisualization(edgePoints, nPoints);
	#endif
	
}

void CResolutionProject::AutoEdgeDetection(const PARE &point, CEdgeTargetNew* target, int iSearchRadius /*радиус*/) const
{
	//length - длина в пикселах 
	//получаем приблизительное (с точностью в 30 градусов положение резкого края) в промежутке [0, PI]
	if(!GetRaster())
		return;
	
	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;

	double dRadius = widthPixel * iSearchRadius;
		
	CIC_Rect3DD geoRect;
	geoRect.top = point.y + iSearchRadius * heightPixel;
	geoRect.bottom = point.y - iSearchRadius * heightPixel;
	geoRect.left = point.x - iSearchRadius * widthPixel;
	geoRect.right = point.x + iSearchRadius * widthPixel;
	
	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;
	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y;

	INT64 iImageHeight = GetRaster()->GetHeight();

	IC_RECT64 boundRect;

	boundRect.left = (_int64)((point.x  - upperLeft_x) / widthPixel - iSearchRadius);
	boundRect.right = (_int64)((point.x  - upperLeft_x) / widthPixel + iSearchRadius);
	boundRect.top = (_int64)((upperLeft_y - point.y) / heightPixel - iSearchRadius);
	boundRect.bottom = (_int64)((upperLeft_y - point.y) / heightPixel + iSearchRadius);

	int nPoints;
	PARE *edgePoints;
	//получаем приближенное значение угла
	double angle = ApproximateAngleDetection(point, (double)iSearchRadius);

	//заполнем массив точек с максимальным градиентом edgePoints
	if( angle < M_PI / 4 || angle > (M_PI * 3) / 4 )
	{
		nPoints = (int)(fabsl((long double)(boundRect.right - boundRect.left))) + 1;
		edgePoints = new PARE[nPoints];
		GetLinesMaxGradPoints(geoRect, edgePoints, nPoints, ROWS);
		GetEdgePosition(target, edgePoints, nPoints, ROWS, m_projectParams.dDefaultTargetLength/* * widthPixel*/, m_projectParams.dDefaultTargetWidth/* * widthPixel*/ /*20 * m_dLocalPixelSizeX*/);
	}
	else
	{
		nPoints = (int)(fabsl((long double)(boundRect.top - boundRect.bottom))) + 1;
		edgePoints = new PARE[nPoints];
		GetLinesMaxGradPoints(geoRect, edgePoints, nPoints, LINES);
		GetEdgePosition(target, edgePoints, nPoints, LINES, m_projectParams.dDefaultTargetLength/* * widthPixel*/, m_projectParams.dDefaultTargetWidth/* * widthPixel*/ /*20 * m_dLocalPixelSizeX*/);
	}

	#ifdef _DEBUG
		//PointsVisualization(edgePoints, nPoints);
	#endif
	
	delete[] edgePoints;
}

void CResolutionProject::GetLinesMaxGradPoints(const CIC_Rect3DD &geoRect, PARE *points, int nPoints, bool bDirection, bool bAbsoluteMaxCalculation) const
{
	// проходим по строкам изображения в прямоугольнике с геоокоординатами, определяемыми geoRect
	// возвращаем массив геокоординат точек с максимальным градиентом
	
	if(!GetRaster())
		return;

	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;

	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y; 
	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;

	IC_RECT64 boundRect;

	//вычисляем координаты описывающего прямоугольника
	boundRect.left = (_int64)((geoRect.left - upperLeft_x) / widthPixel);
	boundRect.right = (_int64)((geoRect.right - upperLeft_x) / widthPixel);
	boundRect.top = (_int64)((upperLeft_y - geoRect.top) / heightPixel);
	boundRect.bottom = (_int64)((upperLeft_y - geoRect.bottom) / heightPixel);


	ASSERT( (bDirection == LINES && nPoints == (int)fabsl((long double)(boundRect.top - boundRect.bottom)) + 1) ||
		    (bDirection == ROWS  && nPoints == (int)fabsl((long double)(boundRect.right - boundRect.left)) + 1) );

	///////////////////////////////////////////////////////////////////
	// Считываем из растра необходимые пикселы 
	///////////////////////////////////////////////////////////////////

	LPBYTE pBlock;
	INT64 iImageWidth, iImageHeight, wEnd, hEnd;
	int iBlockWidth, iBlockHeight, hb, wb, h, w;
	DWORD dwOffsetInBlock;
	
	//DWORD dwPixelSize, dwBlockSize, dwOffsetInBlock;
	GETVALUEPROC pGetValueProc;
	SETVALUEPROC pSetValueProc;
	double dValue,  dMin,  dMax;
	
	CIC_Rect64 rcBlocks;
	CIC_Rect64 rcImage;
	CIC_Rect3DD rcGeo;

	iImageWidth = GetRaster()->GetWidth();
	iImageHeight = GetRaster()->GetHeight();
	iBlockWidth = GetRaster()->GetBlockWidth();
	iBlockHeight = GetRaster()->GetBlockHeight();
	pGetValueProc = CRasterPROCS::PointerGetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));
	pSetValueProc = CRasterPROCS::PointerSetValueProc(GetRaster()->GetPixelType(m_projectParams.iBand));

	GetRaster()->GetImageInfo()->GetDefaultRange(GetRaster()->GetPixelType(m_projectParams.iBand), &dMin, &dMax);
	
	//создаем вспомогательный массив для хранения текущих максимальных значений градиента для каждой строки/столбца
	double *pMaxGradient = new double[nPoints];
	//для подсчета градиента по формуле I[n] - I[n-1] храним массив значений яркостей для предыдущих пикселов каждой строки/столбца
	double *pPreviousValue = new double[nPoints];

	for(int i = 0; i < nPoints; i++)
	{
		pMaxGradient[i] = 0;
		pPreviousValue[i] = 0;
	}

	//////////////////////////////////////////////////////////

	rcBlocks.top	= boundRect.top / iBlockHeight;
	rcBlocks.bottom = boundRect.bottom / iBlockHeight;
	rcBlocks.left	= boundRect.left / iBlockWidth;
	rcBlocks.right	= boundRect.right / iBlockWidth;

	Cdfi_RasterStream rasterStream(GetRaster(), 1);
	INT64 iCurrentHeightPosition, iCurrentWidthPosition;
	double dCurrentGeoHeightPosition, dCurrentGeoWidthPosition;
	
	for(hb = (int)rcBlocks.top; hb <= (int)rcBlocks.bottom; hb++)
	{
		for(wb = (int)rcBlocks.left; wb <= (int)rcBlocks.right; wb++) 
		{
			pBlock = (LPBYTE)rasterStream.Get_lpBlock(0, m_projectParams.iBand, 0, wb, hb);
			ASSERT(pBlock);
			
			wEnd = (wb == rcBlocks.right) ? (boundRect.right - wb*iBlockWidth) : (iBlockWidth - 1);
			hEnd = (hb == rcBlocks.bottom) ? (boundRect.bottom - hb*iBlockHeight) : (iBlockHeight - 1);
			
			for(h = 0; h <= hEnd; h++)
			{
				iCurrentHeightPosition = hb * iBlockHeight + h;
				dCurrentGeoHeightPosition = lowerRight_y + (iImageHeight - 1 - iCurrentHeightPosition) * heightPixel;

				for(w = 0, dwOffsetInBlock = (DWORD)(h*iBlockWidth); w <= wEnd; w++, dwOffsetInBlock++)
				{
					iCurrentWidthPosition = wb * iBlockWidth + w;
					dCurrentGeoWidthPosition = upperLeft_x + iCurrentWidthPosition * widthPixel;
							
					if(iCurrentHeightPosition >= boundRect.top && iCurrentHeightPosition <= boundRect.bottom &&
					   iCurrentWidthPosition >= boundRect.left && iCurrentWidthPosition <= boundRect.right)
					{
						dValue = pGetValueProc((BYTE*)pBlock, dwOffsetInBlock);
						
						if(bDirection == LINES)
						{
							INT64 iCurrentLineIndex = iCurrentHeightPosition - boundRect.top;
							if(iCurrentWidthPosition != boundRect.left)
							{
								double dTmpGradValue = 0;
								//if(bAbsoluteMaxCalculation)
									dTmpGradValue = fabsl(pPreviousValue[iCurrentLineIndex] - dValue);
								//else
								//	dTmpGradValue = pPreviousValue[iCurrentLineIndex] - dValue;

								if(dTmpGradValue > pMaxGradient[iCurrentLineIndex])
								{
									pMaxGradient[iCurrentLineIndex] = dTmpGradValue;
									points[iCurrentLineIndex].y = (double) iCurrentHeightPosition; //dCurrentGeoHeightPosition;
									points[iCurrentLineIndex].x = (double) iCurrentWidthPosition - 0.5; //dCurrentGeoWidthPosition - widthPixel / 2;
								}
							}	
							pPreviousValue[iCurrentLineIndex] = dValue;
						}
						else
						{
							INT64 iCurrentRowIndex = iCurrentWidthPosition - boundRect.left;
							if(iCurrentHeightPosition != boundRect.top)
							{
								//if(bAbsoluteMaxCalculation)
								//	dTmpGradValue = fabsl(pPreviousValue[iCurrentLineIndex] - dValue);
								//else
									double dTmpGradValue = fabsl(pPreviousValue[iCurrentRowIndex] - dValue);
								if(dTmpGradValue > pMaxGradient[iCurrentRowIndex])
								{
									pMaxGradient[iCurrentRowIndex] = dTmpGradValue;
									points[iCurrentRowIndex].y = (double) iCurrentHeightPosition - 0.5; //dCurrentGeoHeightPosition + heightPixel / 2;
									points[iCurrentRowIndex].x = (double) iCurrentWidthPosition; //dCurrentGeoWidthPosition;
								}
							}
							pPreviousValue[iCurrentRowIndex] = dValue;
						}
					}
				}
			}
			GetRaster()->ReleaseLPBlock(m_projectParams.iBand, 0, wb, hb);
		}
	}

	delete[] pMaxGradient;
	delete[] pPreviousValue;
}

void CResolutionProject::GetEdgePosition(CEdgeTargetNew *edgeTarget, const PARE *points, int nPoints, bool bDirection, double dLength, double dWidth) const
{
	if(!GetRaster())
		return;

	double A, B; //коэффициенты прямой, проходящей через резкий край
	
	//Y = A * X + B - метод наименьших квадратов
	
	//        N * SUM(x_i * y_i) - SUM(x_i) * SUM(y_i)              SUM(y_i) - A * SUM(x_i)
	//  A = ----------------------------------------------     B = ------------------------
	//             N * SUM(x_i ^ 2) - SUM(x_i) ^ 2                             N
	
	double dSumX = 0, dSumY = 0, dSumX2 = 0, dSumXY = 0;
	double x, y;

	if(nPoints == 0)
		return;

	double * pX = new double[nPoints];
	double * pY = new double[nPoints];

	PARE centerPoint;
	centerPoint.x = 0;
	centerPoint.y = 0;
  
	for(int i = 0; i < nPoints; i++)
	{
		centerPoint.x += points[i].x;
		centerPoint.y += points[i].y;
	}
	centerPoint.x /= nPoints;
	centerPoint.y /= nPoints;

	PARE pointOnLine; //точка на прямой
	CVector lineDirectionVector;

	if(bDirection == ROWS)
	{
		for(int i = 0; i < nPoints; i++)
		{
			pX[i] = points[i].x - points[0].x;
			pY[i] = points[i].y - points[0].y;
		}

		for(int i = 0; i < nPoints; i++)
		{
			x = pX[i];//points[i].x;
			y = pY[i];//points[i].y;
			dSumX += x;
			dSumY += y;
			dSumX2 += x * x;
			dSumXY += x * y;
		}

		A = (nPoints * dSumXY - dSumX * dSumY) / (nPoints * dSumX2 - dSumX * dSumX);
		B = (dSumY - A * dSumX) / nPoints;

		// вычисляем среднюю точку по всем точкам массива points 
	
		//вычисляем проекцию точки centerPoint на прямую A * X + B

		pointOnLine.x = 0 + points[0].x;
		pointOnLine.y = B + points[0].y;

		//направляющий вектор прямой Ax + B единичной длины
		lineDirectionVector.InitializeXY(1.0, A);
	}
	else if(bDirection == LINES)
	{
		for(int i = 0; i < nPoints; i++)
		{
			pY[i] = points[i].x - points[0].x;
			pX[i] = points[i].y - points[0].y;
		}

		for(int i = 0; i < nPoints; i++)
		{
			x = pX[i];//points[i].x;
			y = pY[i];//points[i].y;
			dSumX += x;
			dSumY += y;
			dSumX2 += x * x;
			dSumXY += x * y;
		}

		A = (nPoints * dSumXY - dSumX * dSumY) / (nPoints * dSumX2 - dSumX * dSumX);
		B = (dSumY - A * dSumX) / nPoints;

		// вычисляем среднюю точку по всем точкам массива points 
		
		pointOnLine.x = B + points[0].x;
		pointOnLine.y = 0 + points[0].y;

		//направляющий вектор прямой Ax + B единичной длины
		lineDirectionVector.InitializeXY(A, 1.0);
	}

	lineDirectionVector /= lineDirectionVector.GetLength();
	//                                  ----------------------->
	//вычисляем длину проекции вектора (pointOnLine, centerPoint) на прямую
	double dProjectionLength = Scalar(CVector(pointOnLine, centerPoint), lineDirectionVector);
	
	// получаем центр тест-объекта
	centerPoint = lineDirectionVector * dProjectionLength + pointOnLine;

	PARE edgePoints[2];

	//получаем точки в координатах снимка
	edgePoints[0] = lineDirectionVector * dLength / 2 + centerPoint; 
	edgePoints[1] = - lineDirectionVector * dLength / 2 + centerPoint; 
	
	//Пересчитываем координаты точек в геокоординаты:
	
	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;

	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y; 

	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;
		
	INT64 iImageWidth = GetRaster()->GetWidth();
	INT64 iImageHeight = GetRaster()->GetHeight();

	edgePoints[0].x = upperLeft_x + edgePoints[0].x * widthPixel;
	edgePoints[0].y = lowerRight_y + (iImageHeight - 1 - edgePoints[0].y) * heightPixel;
	edgePoints[1].x = upperLeft_x + edgePoints[1].x * widthPixel;
	edgePoints[1].y = lowerRight_y + (iImageHeight - 1 - edgePoints[1].y) * heightPixel;

	//PointsVisualization(edgePoints, 2);

	edgeTarget->Initialize(edgePoints[0], edgePoints[1], dWidth / 2, dWidth / 2, GetLocalPixelSizeX(), GetLocalPixelSizeY());

	delete[] pX;
	delete[] pY;
}

void CResolutionProject::GetLinePosition(CLineTarget *edgeTarget, const PARE *points, int nPoints, bool bDirection, double dLength, double dWidth) const
{
	if(!GetRaster())
		return;

	double A, B; //коэффициенты прямой, проходящей через резкий край
	
	//Y = A * X + B - метод наименьших квадратов
	
	//        N * SUM(x_i * y_i) - SUM(x_i) * SUM(y_i)              SUM(y_i) - A * SUM(x_i)
	//  A = ----------------------------------------------     B = ------------------------
	//             N * SUM(x_i ^ 2) - SUM(x_i) ^ 2                             N
	
	double dSumX = 0, dSumY = 0, dSumX2 = 0, dSumXY = 0;
	double x, y;

	if(nPoints == 0)
		return;

	double * pX = new double[nPoints];
	double * pY = new double[nPoints];

	PARE centerPoint;
	centerPoint.x = 0;
	centerPoint.y = 0;
  
	for(int i = 0; i < nPoints; i++)
	{
		centerPoint.x += points[i].x;
		centerPoint.y += points[i].y;
	}
	centerPoint.x /= nPoints;
	centerPoint.y /= nPoints;

	PARE pointOnLine; //точка на прямой
	CVector lineDirectionVector;

	if(bDirection == ROWS)
	{
		for(int i = 0; i < nPoints; i++)
		{
			pX[i] = points[i].x - points[0].x;
			pY[i] = points[i].y - points[0].y;
		}

		for(int i = 0; i < nPoints; i++)
		{
			x = pX[i];//points[i].x;
			y = pY[i];//points[i].y;
			dSumX += x;
			dSumY += y;
			dSumX2 += x * x;
			dSumXY += x * y;
		}

		A = (nPoints * dSumXY - dSumX * dSumY) / (nPoints * dSumX2 - dSumX * dSumX);
		B = (dSumY - A * dSumX) / nPoints;

		// вычисляем среднюю точку по всем точкам массива points 
	
		//вычисляем проекцию точки centerPoint на прямую A * X + B

		pointOnLine.x = 0 + points[0].x;
		pointOnLine.y = B + points[0].y;

		//направляющий вектор прямой Ax + B единичной длины
		lineDirectionVector.InitializeXY(1.0, A);
	}
	else if(bDirection == LINES)
	{
		for(int i = 0; i < nPoints; i++)
		{
			pY[i] = points[i].x - points[0].x;
			pX[i] = points[i].y - points[0].y;
		}

		for(int i = 0; i < nPoints; i++)
		{
			x = pX[i];//points[i].x;
			y = pY[i];//points[i].y;
			dSumX += x;
			dSumY += y;
			dSumX2 += x * x;
			dSumXY += x * y;
		}

		A = (nPoints * dSumXY - dSumX * dSumY) / (nPoints * dSumX2 - dSumX * dSumX);
		B = (dSumY - A * dSumX) / nPoints;

		// вычисляем среднюю точку по всем точкам массива points 
		
		pointOnLine.x = B + points[0].x;
		pointOnLine.y = 0 + points[0].y;

		//направляющий вектор прямой Ax + B единичной длины
		lineDirectionVector.InitializeXY(A, 1.0);
	}

	lineDirectionVector /= lineDirectionVector.GetLength();
	//                                  ----------------------->
	//вычисляем длину проекции вектора (pointOnLine, centerPoint) на прямую
	double dProjectionLength = Scalar(CVector(pointOnLine, centerPoint), lineDirectionVector);
	
	// получаем центр тест-объекта
	centerPoint = lineDirectionVector * dProjectionLength + pointOnLine;

	PARE edgePoints[2];

	//получаем точки в координатах снимка
	edgePoints[0] = lineDirectionVector * dLength / 2 + centerPoint; 
	edgePoints[1] = - lineDirectionVector * dLength / 2 + centerPoint; 
	
	//Пересчитываем координаты точек в геокоординаты:
	
	double upperLeft_x = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.x;
	double upperLeft_y = GetRaster()->GetGeoInfo()->GetMapinfo()->upperLeftCenter.y;

	double lowerRight_y = GetRaster()->GetGeoInfo()->GetMapinfo()->lowerRightCenter.y; 

	double widthPixel  = GetRaster()->GetGeoInfo()->GetPixelSize().width;
	double heightPixel = GetRaster()->GetGeoInfo()->GetPixelSize().height;
		
	INT64 iImageWidth = GetRaster()->GetWidth();
	INT64 iImageHeight = GetRaster()->GetHeight();

	edgePoints[0].x = upperLeft_x + edgePoints[0].x * widthPixel;
	edgePoints[0].y = lowerRight_y + (iImageHeight - 1 - edgePoints[0].y) * heightPixel;
	edgePoints[1].x = upperLeft_x + edgePoints[1].x * widthPixel;
	edgePoints[1].y = lowerRight_y + (iImageHeight - 1 - edgePoints[1].y) * heightPixel;

	//PointsVisualization(edgePoints, 2);

	edgeTarget->Initialize(edgePoints[0], edgePoints[1], dWidth / 2, dWidth / 2, GetLocalPixelSizeX(), GetLocalPixelSizeY());

	delete[] pX;
	delete[] pY;
}


void CResolutionProject::PointsVisualization(const PARE *points, int nPoints) const
{
	Cdfi_Vector *pVector;
	Cdfi_VectorObject *pVObject;
	Cdfi_VectorObjectStream *pVOStream;
	//CIC_StyleObject *style;

	//style->

	IC_POINT3DD *points3DD = new IC_POINT3DD[nPoints];

	for(int i = 0; i < nPoints; i++)
	{
		points3DD[i].x = points[i].x;
		points3DD[i].y = points[i].y;
		points3DD[i].z = 0;
	}

	int iLayer = 0;

	if(GetDocument())
	{
		CIC_VectorInfo vectorInfo(GetCdfi()->GetApp());
		vectorInfo.InitGeoInfo(GetRaster()->GetGeoInfo());

		iLayer = GetDocument()->GetLayers()->AddLayer(&vectorInfo, TRUE);
		
		if(iLayer >= 0)
		{
			pVector = GetDocument()->GetLayers()->GetLayer(iLayer)->GetVector();
			pVObject = pVector->CreateObject(VE_POINT);
			//pVObject->SetStyle( );
			pVObject->GetObjPoint()->AddPoints(nPoints, points3DD);

			pVOStream  = new Cdfi_VectorObjectStream(pVector, VMEM_READWRITE);
			pVOStream->AddObject(pVObject);
		}
	}
	
	pVOStream->Stop();
	if (pVOStream)
		delete pVOStream;
	pVOStream = NULL;
	if (pVObject)
		pVObject->DeleteObject();

	delete[] points3DD;
}

void CResolutionProject::GetInitialESFParameters(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param/*, vector<INTERPOLATION_POINT> &outputPoints*/) const
{
	param.x0 = 0;
	param.N = 2.0; 
	param.R = /*GetLocalPixelSizeX() **/ 2;

	if(points.size() != 0)
	{
		param.minValue = points[0].value;
		param.maxValue = points[0].value;
		for(size_t i = 0; i < points.size(); i++ )
		{
			if(points[i].value < param.minValue )
				param.minValue = points[i].value;
			if(points[i].value > param.maxValue )
				param.maxValue = points[i].value;
		}
	}
	else
	{
		return;
	}

	
	getMinMaxX(points, param.x_min, param.x_max);
	
	AdjustMinMaxValue(points, param);

	GetInitialRnX0Value(points, param.R, param.x0, param);

	AdjustMinMaxValue(points, param);
}

void CResolutionProject::InterpolatePointsLSF(vector<INTERPOLATION_POINT> &points) const
{
	int nMinPoints = 40;
	int nWindow = points.size() / nMinPoints;
	if(nWindow <= 1)
	{
		sort(points.begin(), points.end());
		return;
	}
	else
	{
		sort(points.begin(), points.end());
		int nPointsNew = points.size() / nWindow;
		for(int i = 0; i < nPointsNew; i++)
		{
			double dValue = 0;
			double dX = 0;
			int nPoints = 0;
			for(int j = i * nWindow; j < (i+1) * nWindow; j++)
			{
				dValue += points[j].value;
				dX += points[j].x;
				nPoints++;
			}
			points[i].value = dValue / nPoints;
			points[i].x = dX / nPoints;
		}

		int nTotalPoints = points.size();
		for(int i = 0; i < nTotalPoints -  nTotalPoints / nWindow; i++)
		{
			points.pop_back();
		}
	}
}

void CResolutionProject::InterpolatePoints(vector<INTERPOLATION_POINT> &points) const
{
	double dStepX;
	double dXmin, dXmax;
	getMinMaxX(points, dXmin, dXmax);

	bool bSuccess = false;
	int N = int((dXmax - dXmin) / 0.5) + 10 ; //количество точек дискретизации ESF
	double *pESF = new double[N];
	int *pESFnum = new int[N];

	while(!bSuccess)
	{
		if(N < (dXmax - dXmin) + 10)
		{
			int gugu = 1;
		}

		dStepX = (dXmax - dXmin) / N; //шаг дискретизации

		for(int i = 0; i < N; i++)
		{
			pESF[i] = 0;  
			pESFnum[i] = 0;
		}

		int iSize = points.size();
		for(int i = 0; i < iSize; i++)
		{
			int index = (int)((points[i].x - dXmin) / dStepX);
			index = min(index, N - 1);
			pESF[index] += points[i].value;
			pESFnum[index]++;
			if(index - 1 >= 0)
			{
				pESF[index - 1] += points[i].value;
				pESFnum[index - 1]++;
			}

			if(index + 1 < N)
			{
				pESF[index + 1] += points[i].value;
				pESFnum[index + 1]++;
			}
		}

		bSuccess = true;

		for(int i = 0; i < N; i++)
		{
			if(pESFnum[i] != 0)
				pESF[i] /= pESFnum[i];
			else
			{
				N = (int) ((double)(N) / 1.05);
				bSuccess = false;
				break;
			}
		}
	}

	INTERPOLATION_POINT tmpPoint;
	points.clear();
	
	for(int i = 0; i < N; i++)
	{
		tmpPoint.x = dXmin + (i + 1) * dStepX;
		tmpPoint.value = pESF[i];

		points.push_back(tmpPoint);
	}

	delete[] pESFnum;
	delete[] pESF;

}

void CResolutionProject::GetInitialLSFParameters(vector<INTERPOLATION_POINT> &points, ESF_PARAMETERS_NEW& param/*, vector<INTERPOLATION_POINT> &outputPoints*/) const
{
	getMinMaxX(points, param.x_min, param.x_max);
	
	double dStepX = 0;

	/*ofstream resultingFile;
	resultingFile.open("C:\\temp\\POINTS3.txt", std::ios::out);
	OutputListOfPoints(points, resultingFile);
	resultingFile.close(); 
	*/
	
	param.x0 = 0;
	param.N = 2.0; 
	param.R = 2;

	if(points.size() != 0)
	{
		param.minValue = points[0].value;
		param.maxValue = points[0].value;
		for(size_t i = 0; i < points.size(); i++ )
		{
			if(points[i].value < param.minValue )
				param.minValue = points[i].value;
			if(points[i].value > param.maxValue )
				param.maxValue = points[i].value;
		}
	}
	else
	{
		return;
	}

	
	
	/*
	AdjustMinMaxValue(points, param);

	GetInitialRnX0Value(points, param.R, param.x0, param);

	AdjustMinMaxValue(points, param);
	*/
}



void CResolutionProject::GetInitialRnX0Value(vector<INTERPOLATION_POINT> &points, double &R, double &x0, const ESF_PARAMETERS_NEW& param/*, vector<INTERPOLATION_POINT> &outputPoints*/) const
{
	//предполагается, что x_min, x_max, minValue, maxValue - уже ЗАПОЛНЕНЫ, а param.x0 = 0
	/*	
	int i;
	double dStepX;
	bool bSuccess = false;
	int N = int((param.x_max - param.x_min) / 0.5) + 10 ; //количество точек дискретизации ESF
	double *pESF = new double[N];
	int *pESFnum = new int[N];

	while(!bSuccess)
	{
		dStepX = (param.x_max - param.x_min) / N; //шаг дискретизации

		for(i = 0; i < N; i++)
		{
			pESF[i] = 0;
			pESFnum[i] = 0;
		}

		int iSize = points.size();
		for(i = 0; i < iSize; i++)
		{
			int index = (int)((points[i].x - param.x_min) / dStepX);
			index = min(index, N - 1);
			pESF[index] += points[i].value;
			pESFnum[index]++;
		}

		bSuccess = true;

		for(i = 0; i < N; i++)
		{
			if(pESFnum[i] != 0)
				pESF[i] /= pESFnum[i];
			else
			{
				N = (int) ((double)(N) / 1.05);
				bSuccess = false;
				break;
			}
		}
	}
	
	*/
	double x_09 = 0, x_01 = 0;

	double dStepX = (param.x_max - param.x_min) / points.size(); //шаг дискретизации

	for(int i = 0; i < (int)points.size() ; i++)
	{
	//for(i = (int)((- param.x_min) / dStepX); i < N ; i++)
		if(points[i].value > param.minValue + 0.9 * (param.maxValue - param.minValue))
		{
			x_09 = param.x_min + dStepX / 2 + i * dStepX;
			break;
		}
	}
	 
	for(int i = points.size() - 1; i > 0 ; i--)
	{
	//for(i = (int)((- param.x_min) / dStepX); i > 0  ; i--)
		if(points[i].value < param.minValue + 0.1 * (param.maxValue - param.minValue))
		{
			x_01 = param.x_min - dStepX / 2 + i * dStepX;
			break;
		}
	}
	
	R = abs(x_09 - x_01);

	if(R < dStepX)
		R = dStepX;

	x0 = (x_09 + x_01) / 2;
}


double CResolutionProject::CSFValue(double dFrequency_meters) const
{
	////////////// dFrequency_meters - частота в 1/м

	double dFrequency = dFrequency_meters; // * GetLocalPixelSizeX();

	double result = 0;

	//double C = 0.1 * 255;
	CSF_PARAMETERS_SNR params = m_projectParams.CSFparamsSNR;
	//result = C * (params.dSNRLimit * params.dNoise * dFrequency) / ((params.dMax - params.dMin) * m_projectParams.dContrastLimit * sqrt(params.nTestLines * 2 - 1));
	if(	m_projectParams.CSFparamsSNR.dNoise < std::numeric_limits<double>::epsilon() ||
		params.dNoise < std::numeric_limits<double>::epsilon() ||
		params.dNoise < std::numeric_limits<double>::epsilon() ||
		(1 - m_projectParams.dContrastLimit) < std::numeric_limits<double>::epsilon() ||
		m_projectParams.dMeanPixelValue < std::numeric_limits<double>::epsilon()) 
	{
		return -1;
	}

	result = 4.0 * (dFrequency * M_PI * M_PI * params.dNoise * params.dSNRLimit) / ( 8 * sqrt(sqrt(params.nTestLines * 2 - 1)) * m_projectParams.dMeanPixelValue * m_projectParams.dContrastLimit);
	
	//2 * (params.dSNRLimit * params.dNoise * dFrequency * M_PI * M_PI) /
		//(4 * m_projectParams.dMeanPixelValue * (2 * m_projectParams.dContrastLimit) / (1 - m_projectParams.dContrastLimit) * sqrt(params.nTestLines * 2 - 1));
	
	double coefficient = (M_PI * M_PI * params.dNoise * params.dSNRLimit) / ( 8 * sqrt(sqrt(params.nTestLines * 2 - 1)) * m_projectParams.dMeanPixelValue * m_projectParams.dContrastLimit);
	double test = params.dNoise / m_projectParams.dMeanPixelValue;

	return result;
}

BOOL CResolutionProject::Draw(Cdfi_View *pView, Cdfi_MemoryDC *pMemDC)
{
	int nEdges = m_listOfTargets.size();
	for(int i = 0; i < nEdges; i++)
		m_listOfTargets[i]->Draw(pView, pMemDC);

	return TRUE;
}

int CResolutionProject::GetBoundRect(CIC_Rect3DD& rect)
{
	CIC_PolygonD polygon;
	int nTargets = GetNumberOfTargets();
	if(nTargets != 0)
	{
		GetTarget(0)->GetPolygon(polygon);
		polygon.UpdateBoundRectGeo();
		rect.bottom = polygon.GetBoundRect()->bottom;
		rect.top = polygon.GetBoundRect()->top;
		rect.left = polygon.GetBoundRect()->left;
		rect.right = polygon.GetBoundRect()->right;
		
		for (int i = 0; i < nTargets; i++)
		{
			GetTarget(i)->GetPolygon(polygon);
			polygon.UpdateBoundRectGeo();
			rect.bottom = min(rect.bottom, polygon.GetBoundRect()->bottom);
			rect.top = max(rect.top, polygon.GetBoundRect()->top);
			rect.left = min(rect.left, polygon.GetBoundRect()->left);
			rect.right = max(rect.right, polygon.GetBoundRect()->right);
		}
		return 1;
	}
	else
		return -1;
}

int CResolutionProject::GetBoundRectSelected(CIC_Rect3DD& rect)
{
	CIC_PolygonD polygon;
	int nTargets = GetNumberOfTargets();
	if(nTargets != 0)
	{
		int result = -1;
		
		for (int i = 0; i < nTargets; i++)
		{
			if(GetTarget(i)->GetSelectionStatus())
			{
				if(result == -1)
				{
					GetTarget(i)->GetPolygon(polygon);
					polygon.UpdateBoundRectGeo();
					rect.bottom = polygon.GetBoundRect()->bottom;
					rect.top = polygon.GetBoundRect()->top;
					rect.left = polygon.GetBoundRect()->left;
					rect.right = polygon.GetBoundRect()->right;
				}
				else
				{
					GetTarget(i)->GetPolygon(polygon);
					polygon.UpdateBoundRectGeo();
					rect.bottom = min(rect.bottom, polygon.GetBoundRect()->bottom);
					rect.top = max(rect.top, polygon.GetBoundRect()->top);
					rect.left = min(rect.left, polygon.GetBoundRect()->left);
					rect.right = max(rect.right, polygon.GetBoundRect()->right);
				}

				result = 1;
			}
		}
		return result;
	}
	else
		return -1;
}

double CResolutionProject::GetContrastLimit() const
{
	return m_projectParams.dContrastLimit;
}

double CResolutionProject::GetPixelSizeX() const
{
	return m_projectParams.dPixelSizeX;
}

double CResolutionProject::GetPixelSizeY() const
{
	return m_projectParams.dPixelSizeY;
}

double CResolutionProject::GetDefaultTargetLength() const
{
	return m_projectParams.dDefaultTargetLength;
}

double CResolutionProject::GetDefaultTargetWidth() const
{
	return m_projectParams.dDefaultTargetWidth;
}

void CResolutionProject::SetDefaultTargetLength(double length) 
{
	m_projectParams.dDefaultTargetLength = length;
}

void CResolutionProject::SetDefaultTargetWidth(double width)
{
	m_projectParams.dDefaultTargetWidth = width;
}

void CResolutionProject::SetLocalPixelSizeX(double sizeX)
{
	m_projectParams.dLocalPixelSizeX = sizeX;
}

void CResolutionProject::SetLocalPixelSizeY(double sizeY)
{
	m_projectParams.dLocalPixelSizeY = sizeY;
}


double CResolutionProject::GetLocalPixelSizeX() const
{
	if(GetRaster())
	{
		double dLocalPixelSizeX = GetRaster()->GetGeoInfo()->GetPixelSize().width; // размер пиксела в локальных единицах измерения
		//m_projectParams.dLocalPixelSize = dLocalPixelSizeX;
		return dLocalPixelSizeX;
	}
	else
	{	
		return m_projectParams.dLocalPixelSizeX;
	}
}
double CResolutionProject::GetLocalPixelSizeY() const
{
	if(GetRaster())
	{
		double dLocalPixelSizeY = GetRaster()->GetGeoInfo()->GetPixelSize().height;
		return dLocalPixelSizeY;
	}
	else
	{
		return m_projectParams.dLocalPixelSizeY;
	}
}

int CResolutionProject::HitTest(CPare& ptGeo, int& hitEdgeIndex)
{
	IC_POINT3DD point;
	
	point.x = ptGeo.x;
	point.y = ptGeo.y;
	point.z = 0;
	
	int nEdges = m_listOfTargets.size();
	for(int i = 0; i < nEdges; i++)
	{
		if(m_listOfTargets[i]->HitTest(&point/*, GetLocalPixelSizeX()*/) == INSIDE_TARGET)
		{
			hitEdgeIndex = i; 
			return INSIDE_TARGET;
		}
	}

	hitEdgeIndex = -1;
	return NO_TARGET;	
}

long CResolutionProject::Serialize(void* pMemory, BOOL IsReading)
{
	BYTE *p = (BYTE*)pMemory;
	long size = 0;
	int i;

	if ( !pMemory )
	{
		size += sizeof(m_projectParams);

		size += sizeof(int);

		for(int k = 0; k < (int)(m_listOfTargets.size()); k++)
		{
			size += m_listOfTargets[k]->Serialize(NULL, IsReading);
		}

		size += sizeof(int);

		for(int j = 0; j < (int)(m_vectorNoise.size()); j++)
		{
			size += sizeof(m_vectorNoise[j]);
		}

		return size;
	}
	else
	{
		if (IsReading)
		{
			memcpy( &m_projectParams, p, sizeof(m_projectParams) );
			p += sizeof(m_projectParams);

			int N;

			memcpy( &N, p, sizeof(N) );
			p += sizeof(N);
			m_listOfTargets.resize(N);
			
			for ( i = 0; i < N; i++ )
			{
				int iType;
				memcpy( &iType, p, sizeof(int));

				if(iType == LINE_TARGET)
				{
					m_listOfTargets[i] = new CLineTarget();
					p += m_listOfTargets[i]->Serialize(p, IsReading);
				}
				else if(iType == EDGE_TARGET)
				{
					m_listOfTargets[i] = new CEdgeTargetNew();
					p += m_listOfTargets[i]->Serialize(p, IsReading);
				}
				//p += m_listOfTargets[i]->Serialize(p, IsReading);
			}

			memcpy( &N, p, sizeof(N) );
			p += sizeof(N);
			m_vectorNoise.resize(N);
			
			for ( i = 0; i < N; i++ )
			{
				memcpy( &m_vectorNoise[i], p, sizeof(m_vectorNoise[i]) );
				p += sizeof(m_vectorNoise[i]);
			}
		}
		else
		{
			memcpy( p, &m_projectParams, sizeof(m_projectParams) );
			p += sizeof(m_projectParams);
		
			int N = m_listOfTargets.size();
			
			memcpy( p, &N, sizeof(N) );
			p += sizeof(N);
		 
			for ( i = 0; i < N; i++ )
			{
				int size = m_listOfTargets[i]->Serialize(p, IsReading);
				p += size;
			}

			N = m_vectorNoise.size();

			memcpy( p, &N, sizeof(N) );
			p += sizeof(N);

			for ( i = 0; i < N; i++ )
			{
				memcpy( p, &m_vectorNoise[i], sizeof(m_vectorNoise[i]) );
				p += sizeof(m_vectorNoise[i]);
			}
		}
	}

	return long(p - (BYTE*)pMemory);
}

Cdfi_Document *CResolutionProject::GetDocument() const
{
	if(m_pResolutionLayer)
	{
		Cdfi_Document *document = ((Cdfi_LayerObject*) m_pResolutionLayer)->GetDocument();
		return document;
	}
	else
		return NULL;
}

Cdfi_Raster *CResolutionProject::GetRaster() const
{
	Cdfi_Document *pDocument = GetDocument();
	
	int iRasterIndex = m_projectParams.iRasterIndex;

	if(iRasterIndex == -1)
		return NULL;
	
	if(pDocument)
	{
		BOOL bLockedWrite = pDocument->IsLockDocument(MODE_ACCESS_WRITE);
		if(bLockedWrite)
			pDocument->UnlockDocument(MODE_ACCESS_WRITE);
		BOOL bLocked = pDocument->IsLockDocument(MODE_ACCESS_READ);
		if(!bLocked)
			pDocument->LockDocument(MODE_ACCESS_READ);
		
		if(pDocument->IsLockDocument(MODE_ACCESS_READ))
		{
			Cdfi_Raster * raster = pDocument->GetLayers()->GetLayer(m_projectParams.iRasterIndex)->GetRaster();
			
			if(!bLocked && pDocument->IsLockDocument(MODE_ACCESS_READ))
				pDocument->UnlockDocument(MODE_ACCESS_READ);
			if(bLockedWrite && !pDocument->IsLockDocument(MODE_ACCESS_WRITE))
				pDocument->LockDocument(MODE_ACCESS_WRITE);

			return raster;
		}
		else 
		{
			if(!bLocked && pDocument->IsLockDocument(MODE_ACCESS_READ))
				pDocument->UnlockDocument(MODE_ACCESS_READ);
			if(bLockedWrite && !pDocument->IsLockDocument(MODE_ACCESS_WRITE))
				pDocument->LockDocument(MODE_ACCESS_WRITE);
			return NULL;
		}
	}
	else 
		return NULL;
}

CString CResolutionProject::GetDocName() const
{
	Cdfi_Document *pDocument = GetDocument();
	if(pDocument)
	{
		CString title = pDocument->GetTitle();
		return title;
	}
	else
		return _T("");
}

CString CResolutionProject::GetBandName() const
{
	Cdfi_Raster *raster = GetRaster();
	if(raster)
	{
		CString sBandName = raster->GetImageInfo()->GetBandInfo(m_projectParams.iBand)->BandName;
		return sBandName;
	}
	else
		return _T("");
}
CString CResolutionProject::GetRasterLayerName() const
{
	Cdfi_Document *pDocument = GetDocument();
	if(pDocument)
	{
		BOOL bLockedWrite = pDocument->IsLockDocument(MODE_ACCESS_WRITE);
		if(bLockedWrite)
			pDocument->UnlockDocument(MODE_ACCESS_WRITE);
		BOOL bLocked = pDocument->IsLockDocument(MODE_ACCESS_READ);
		if(!bLocked)
			pDocument->LockDocument(MODE_ACCESS_READ);
		
		if(pDocument->IsLockDocument(MODE_ACCESS_READ))
		{
			CString sRasterLayerName;
			if(m_projectParams.iRasterIndex != -1)
				sRasterLayerName = pDocument->GetLayers()->GetLayer(m_projectParams.iRasterIndex)->GetLayerName();
			
			else
				sRasterLayerName = "";

			if(!bLocked && pDocument->IsLockDocument(MODE_ACCESS_READ))
				pDocument->UnlockDocument(MODE_ACCESS_READ);
			if(bLockedWrite && !pDocument->IsLockDocument(MODE_ACCESS_WRITE))
				pDocument->LockDocument(MODE_ACCESS_WRITE);

			return sRasterLayerName;
		}
		else
		{
			if(!bLocked && pDocument->IsLockDocument(MODE_ACCESS_READ))
				pDocument->UnlockDocument(MODE_ACCESS_READ);
			if(bLockedWrite && !pDocument->IsLockDocument(MODE_ACCESS_WRITE))
				pDocument->LockDocument(MODE_ACCESS_WRITE);
			return _T("");
		}
	}
	else
		return _T("");
}

CString CResolutionProject::GetResolutionLayerName() const
{
	Cdfi_Document *pDocument = GetDocument();

	if(m_pResolutionLayer /*&& pDocument*/)
	{
		BOOL bLockedWrite = pDocument->IsLockDocument(MODE_ACCESS_WRITE);
		if(bLockedWrite)
			pDocument->UnlockDocument(MODE_ACCESS_WRITE);
		BOOL bLocked = pDocument->IsLockDocument(MODE_ACCESS_READ);
		if(!bLocked)
			pDocument->LockDocument(MODE_ACCESS_READ);
		
		if(pDocument->IsLockDocument(MODE_ACCESS_READ))
		{
			CString sResolutionLayerName = ((Cdfi_LayerObject*) m_pResolutionLayer)->GetLayer()->GetLayerName();
		
			if(!bLocked && pDocument->IsLockDocument(MODE_ACCESS_READ))
				pDocument->UnlockDocument(MODE_ACCESS_READ);
			if(bLockedWrite && !pDocument->IsLockDocument(MODE_ACCESS_WRITE))
				pDocument->LockDocument(MODE_ACCESS_WRITE);

			return sResolutionLayerName;
		}
		else
		{
		
			return _T("");
		}
	}
	else
		return _T("");
}
 
int CResolutionProject::GetBand() const
{
	return m_projectParams.iBand;
}

int CResolutionProject::GetRasterIndex() const
{
	return m_projectParams.iRasterIndex;
}

void CResolutionProject::SetRasterIndex(int iRasterInd)
{
	m_projectParams.iRasterIndex = iRasterInd;
}

void CResolutionProject::SetBand(int iBand)
{
	m_projectParams.iBand = iBand;
}

CDFI *CResolutionProject::GetCdfi() const
{
	Cdfi_Document *pDocument = GetDocument();
	if(pDocument)
		return pDocument->GetMainFrame()->GetDfi();
	else
		return NULL;
}



PROJECT_PARAMETERS CResolutionProject::GetParams() const
{
	return m_projectParams;
}

void CResolutionProject::SetResolutionLayer(CResolutionLayer *pResolutionLayer)
{
	m_pResolutionLayer = pResolutionLayer;
}

CResolutionLayer *CResolutionProject::GetResolutionLayer()
{
	return m_pResolutionLayer;
}

int CResolutionProject::GetNumberOfNoiseTargets()
{
	return m_vectorNoise.size();
}

CNoiseTarget *CResolutionProject::GetNoiseTarget(int i)
{
	int iSize = m_vectorNoise.size();
	if(i >= iSize)
		return NULL;

	return &m_vectorNoise[i]; 
}