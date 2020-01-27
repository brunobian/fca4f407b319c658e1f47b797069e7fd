modelo <- function(str, v, data, n1){
  print(c(str,v))
  if (n1){
    if (str == 'w2v'){
      data['este'] = data[paste0('w2v',v,'.distancia.conSW_next')]  
    } else if (str == 'FT_w'){
      data['este'] = data[paste0('FT',v,'.distancia.promedio_conSW_wiki_next')]  
    } else if (str == 'FT_c'){
      data['este'] = data[paste0('FT',v,'.distancia.promedio_conSW_corpus_next')]  
    } else{
      data['este'] = data[paste0('LSA',v,'.', str,'.conSW_next')]
    }
    corr <- cor(data["este"], data$CLOZE_pred_next, method = c("spearman"))
  } else {
    if (str == 'w2v'){
      data['este'] = data[paste0('w2v',v,'.distancia.conSW')]  
    } else if (str == 'ngram'){
      data['este'] = data[paste0(v)]  
    } else if (str == 'FT_w'){
      data['este'] = data[paste0('FT',v,'.distancia.promedio_conSW_wiki')]  
    } else if (str == 'FT_c'){
      data['este'] = data[paste0('FT',v,'.distancia.promedio_conSW_corpus')]  
    } else{
      data['este'] = data[paste0('LSA',v,'.', str,'.conSW')]
    }
    corr <- cor(data["este"], data$CLOZE_pred, method = c("spearman"))
  }

  data['este'] <- scale(data['este'], scale=FALSE, center=TRUE)
  # print('escaleado')
  # Corro modelo
  M <- lmer(logFPRT ~ Nlaunchsite + invlength*freq + 
              rpl + rpt + rps + 
              este +
              (1|sujid) + (1|textid) + (1|wordid), data = data, REML = FALSE)
  aic  <- extractAIC(M)[2]
  tval <- coef(summary(M))[,"t value"]['este']
  
  return(c(aic, tval, corr))
}

