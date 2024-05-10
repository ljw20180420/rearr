<script setup>
import uploadFile from './components/uploadFile.vue';
import runJob from './components/runJob.vue';
import setValues from './components/setValues.vue';
</script>

<template>
  <div class="outer">
    <img src="./assets/images/projectLogic.png">
  </div>

  <div class="outer">
    <uploadFile url="targetFile" description="target file" />
    <uploadFile url="pairFile" description="pair file" />
    <runJob url="removeDup" description="remove duplicate file" />
  </div>

  <div class="outer">
    <uploadFile url="targetSpliter" description="target spliter" />
    <uploadFile url="pairSpliter" description="pair spliter" />
    <runJob url="bowtie2build" description="bowtie2 build index" />
  </div>

  <div class="outer">
    <uploadFile url="rmDupFile" description="file without duplicates" />
    <setValues :params='["minScoreTarget", "minScorePair"]' :descriptions='["minimal alignment score of target spliter", "minimal alignment score of pair spliter"]' :values='[null, null]' />
    <runJob url="demultiplex" description="demultiplex" />
  </div>

  <div class="outer">
    <uploadFile url="toAlignFile" description="file to align" />
    <uploadFile url="refFile" description="reference file" />
    <setValues :params='["s0", "s1", "s2", "u", "v", "ru", "rv", "qu", "qv"]' :descriptions='["mismatching score", "matching score for non-extension reference part", "matching score for extension reference part", "gap-extending penalty", "gap-opening penalty", "gap-extending penalty for unaligned reference end", "gap-opening penalty for unaligned reference end", "gap-extending penalty for unaligned query part", "gap-opening penalty for unaligned query part"]' :values='[-6, 4, 2, -3, -9, 0, 0, 0, -5]' :selectParams='["PAM1", "PAM2"]' :selectDescriptions='["PAM for cut1", "PAM for cut2"]' :selectValues='[["NGG", "CCN"], ["NGG", "CCN"]]' />
    <runJob url="rearrange" description="align reads and correct microhomologies" />
  </div>

  <div class="outer">
    <h1>Monitor Progresses</h1>
    <a href="/flower" target="_blank">Celery Flower</a>
  </div>

  <div class="outer">
    <h1>Shiny Applications</h1>
    <a href="/shiny" target="_blank">Shiny Server</a>
  </div>
</template>

<style>
    .outer {
        outline: 0px solid #CCC;
        border: 10px solid #DDD;
        background-color: #999;
        overflow: auto;
    }
</style>
