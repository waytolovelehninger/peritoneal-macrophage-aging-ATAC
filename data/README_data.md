# 데이터 받는 방법

이 저장소에는 데이터 파일이 들어 있지 않습니다. 아래 파일을 직접 내려받아 `data/` 폴더에 넣은 뒤 스크립트를 실행하세요.

## 필요한 파일: `data/media-3.xlsx`

| 항목 | 내용 |
|---|---|
| 출처 | Lu et al., bioRxiv 2025.06.11.659200 (v1, 2025-06-17 게시) |
| 파일 | bioRxiv supplementary의 **`media-3.xlsx`** = **Supplementary Table S2** |
| 내려받는 곳 | 논문 페이지 → *Supplementary Material* 탭: <https://www.biorxiv.org/content/10.1101/2025.06.11.659200v1.supplementary-material> |
| 저장 위치 | `data/media-3.xlsx` (파일 이름 그대로) |
| 크기 | 31,724,858 bytes (약 30 MB) |
| MD5 | `43B7E916BD80874A12A5E42FCAF52487` |
| SHA-256 | `5D400048E3A4063F7ECCDB31F31AEAA2F950B7B6348ABB4F8EAB2E60B6EE7310` |

파일 안의 시트 구성:

| 시트 | 내용 | 이 분석에서 사용 |
|---|---|---|
| Index | 시트 설명 | - |
| A | RNA-seq, 유전자별 DESeq2 결과 | 사용 안 함 |
| **B** | **ATAC-seq, peak별 DESeq2 결과** | **사용** |
| C | H3K4me3, peak별 DESeq2 결과 | 사용 안 함 |

### 파일 이름 주의

bioRxiv의 **`media-2.xlsx`는 Supplementary Table S1**(실험에 쓴 쥐 정보)이고, ATAC-seq 결과가 든 Table S2는 **`media-3.xlsx`**입니다. 이름이 비슷하니 헷갈리지 않도록 주의하세요.

### 내려받은 파일 확인 (선택)

파일이 올바른지 MD5 값으로 확인할 수 있습니다. 위 표의 값과 같아야 합니다.

```sh
# Windows PowerShell
Get-FileHash data/media-3.xlsx -Algorithm MD5

# macOS
md5 data/media-3.xlsx

# Linux
md5sum data/media-3.xlsx
```

## 참고: 쥐 나이 정보 (Supplementary Table S1, `media-2.xlsx`)

스크립트 실행에는 필요 없습니다. README에 적은 Table S1 계산값(young 6.54개월, old 22.57개월)은 이 파일의 sheet A("Aging C57BL/6JNia mice information")에서 `Experiment`가 "Peritoneal macrophage ATAC-seq/purity"인 20마리(성별·나이 그룹당 5마리)의 생년월일(`Date of Birth`)과 안락사일(`Date of Euthanaisa`, 원문 철자 그대로)로 계산했습니다.

- young: 2020-12-15 출생 → 2021-07-02 안락사 (6.54개월)
- old: 2019-08-15 출생 → 2021-07-02 안락사 (22.57개월)
- 한 달 = 30.4375일로 계산

이 값은 논문 본문의 나이 표기(young 4-5개월, old 20-21개월)와 다르며, 그 이유는 **확인 필요**입니다. 두 기준 모두 두 그룹의 차이가 약 16개월이므로 ×16 환산과 p값·padj에는 영향이 없습니다.

## raw data

원본 시퀀싱 데이터는 NCBI BioProject **PRJNA1095555**에 등록되어 있습니다. 이 저장소는 raw data를 다시 분석하지 않고, 원 논문 저자들이 공개한 DESeq2 결과표(Table S2)만 사용합니다.
