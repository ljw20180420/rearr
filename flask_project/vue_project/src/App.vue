<script setup>
import uploadFile from './uploadFile.vue'
</script>

<template>
  <div class="outer">
    <uploadFile url="R1" />
  </div>
  <div class="outer">
    <uploadFile url="removeDup" />
  </div>

  <!-- <div class="outer">
    <h1>Align Reads</h1>
    <AlignReads />
  </div> -->

  <!-- <div class="outer">
    <h1>Monitor Progresses</h1>
    <a href="/flower" target="_blank">Celery Flower</a>
  </div>

  <div class="outer">
    <h1>Shiny Applications</h1>
    <a href="/shiny" target="_blank">Shiny Server</a>
  </div> -->
</template>

<style>
    .outer {
       outline: 0px solid #CCC;
       border: 50px solid #FFF;
       background-color: #999;
    }
</style>
