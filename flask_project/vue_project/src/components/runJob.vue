<script setup>
import { ref } from 'vue';
import axios from 'axios';
const taskId = ref(null);
const props = defineProps({
    'description': String,
    'url': String
});
async function runJob() {
    try{
        const response = await axios.get("/runJob/" + props.url);
        taskId.value = response.data;
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <span>{{ description }}</span><br>
        <button @click="runJob">run</button>
        <a :href="'/download/' + taskId">{{ taskId }}</a>
    </div>
</template>