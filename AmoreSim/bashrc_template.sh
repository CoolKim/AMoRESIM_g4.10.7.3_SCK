#!/bin/bash

# CMAKE를 이용한 CupSoftware Building에 필요한 환경 변수들
# MUON_PATH: Muon flux data가 있는 파일의 경로 (muon simulation을 하지 않는다면 아무 의미없는 문자열을 사용해도 무방)
# NEUTRON_PATH: Neutron flux data가 있는 파일의 경로 (internal simulation을 하지 않는다면 아무 의미없는 문자열을 사용해도 무방)
# SIMOUT_PATH: 시뮬레이션 출력 파일 경로 (잡 클러스터에서 실제 시뮬레이션을 하지 않는다면 아무 의미없는 문자열을 사용해도 무방)
# MUON_CONTOUR_PATH: Muon flux data파일의 생성을 위한 Contour 데이터 파일이 있는 경로 (역시 아무 의미없는 문자열을 사용해도 무방)
# MCOBJS_PATH: MCObjs의 라이브러리 경로 (시뮬레이션 프로그램의 실행을 위해 필요할 수도 있음)

# 환경 변수 설정예시
export MUON_PATH="N/A"
export NEUTRON_PATH="N/A"
export SIMOUT_PATH="N/A"
export MUON_CONTOUR_PATH="N/A"
export MCOBJS_PATH="N/A"
